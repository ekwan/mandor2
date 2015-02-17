import java.util.*;
import java.io.*;
import java.util.concurrent.*;
import java.util.concurrent.atomic.*;
import com.google.common.collect.*;
import com.google.common.util.concurrent.AtomicDouble;
import org.jgrapht.*;
import org.jgrapht.graph.*;
import org.jgrapht.alg.*;
import org.apache.commons.math3.geometry.euclidean.threed.*;

/**
 * Performs various kinds of doubles DEE.
 */
public class DEEDoubles
{
    /** the self-energy of the sidechain and the interaction energy with the backbone, corrected by reference energy */
    public final Map<Rotamer,Double> rotamerSelfEnergies;

    /** the interaction energies between the two rotamers */
    public final Map<RotamerPair,Double> rotamerInteractionEnergies;

    /** the space of all Rotamers allowed at the various positions on the amino acid.  */
    public final List<List<Rotamer>> rotamerSpace;

    /** the incompatible pairs we have at the start of the calculation */
    public final Set<RotamerPair> oldIncompatiblePairs;

    /** any incompatible pairs added during the calculation */
    public final Set<RotamerPair> newIncompatiblePairs;

    /** map from rotamer pairs to epsilon_min */
    public final ConcurrentHashMap<RotamerPair,Double> e_min_map;

    /** map from rotamer pairs to epsilon_max */
    public final ConcurrentHashMap<RotamerPair,Double> e_max_map;

    /** constructor */
    public DEEDoubles(Map<Rotamer,Double> rotamerSelfEnergies,
                      Map<RotamerPair,Double> rotamerInteractionEnergies,
                      List<List<Rotamer>> rotamerSpace,
                      Set<RotamerPair> oldIncompatiblePairs)
    {
        this.rotamerSelfEnergies = rotamerSelfEnergies;
        this.rotamerInteractionEnergies = rotamerInteractionEnergies;
        this.rotamerSpace = rotamerSpace;
        this.oldIncompatiblePairs = oldIncompatiblePairs;
        // effectively a concurrent hash set
        this.newIncompatiblePairs = Collections.newSetFromMap(new ConcurrentHashMap<RotamerPair,Boolean>());
        e_min_map = new ConcurrentHashMap<>();
        e_max_map = new ConcurrentHashMap<>();
    }

    /**
     * Static factory method for creating a DEEDoubles object from a DEESingles object.
     * Rotamer interaction energies between pairs that include a rotamer that is no longer
     * in the rotamer space are removed.  Similarly, single rotamer energies for pruned rotamers are removed.
     */
    public static DEEDoubles createDEEDoubles(DEESingles singles)
    {
        //System.out.print("Creating DEEDoubles objects from DEESingles object...");
        // prune eliminated rotamers from the single rotamer energies
        HashSet<Rotamer> includedRotamers = new HashSet<>();
        for (List<Rotamer> list : singles.rotamerSpace)
            includedRotamers.addAll(list);

        Map<Rotamer,Double> newRotamerSelfEnergies = new HashMap<>(singles.rotamerSelfEnergies);
        newRotamerSelfEnergies.keySet().retainAll(includedRotamers);

        // only retain interaction energies we actually need
        Map<RotamerPair,Double> newRotamerInteractionEnergies = new HashMap<>();
        for (RotamerPair pair : singles.rotamerInteractionEnergies.keySet())
            {
                if ( includedRotamers.contains(pair.rotamer1) && 
                     includedRotamers.contains(pair.rotamer2)    )
                    {
                        Double energy = singles.rotamerInteractionEnergies.get(pair);
                        newRotamerInteractionEnergies.put(pair,energy);
                    }
            }

        //System.out.println("done creating object.");
        return new DEEDoubles(newRotamerSelfEnergies,
                              newRotamerInteractionEnergies,
                              singles.rotamerSpace,
                              singles.incompatiblePairs);
    }

    /**
     * Makes a DEEDoubles object from another one.  Prunes interaction energies for new incompatible pairs (assumes
     * that old incompatible pairs have already been pruned).  Combines incompatible lists.
     * @return new DEEDoubles with 
     */
    public DEEDoubles prune()
    {
        // get rid of incompatible energies from rotamerInteractionEnergies
        // assumes the oldIncompatiblePairs will not be in the keys of rotamerInteractionEnergies
        //System.out.println("Pruning excess pairs from the set of interaction energies...");
        
        int estimatedSize = rotamerInteractionEnergies.size() - newIncompatiblePairs.size();
        Map<RotamerPair,Double> newRotamerInteractionEnergies = new HashMap<>(estimatedSize);
        for (RotamerPair pair : rotamerInteractionEnergies.keySet())
            {
                if ( newIncompatiblePairs.contains(pair) )
                    continue;
                Double energy = rotamerInteractionEnergies.get(pair);
                newRotamerInteractionEnergies.put(pair,energy);
            }

        // combine the incompatible pairs lists
        //System.out.println("Combining incompatible lists...");
        
        estimatedSize = oldIncompatiblePairs.size() + newIncompatiblePairs.size();
        Set<RotamerPair> combinedIncompatiblePairs = new HashSet<>(estimatedSize);
        combinedIncompatiblePairs.addAll(oldIncompatiblePairs);
        combinedIncompatiblePairs.addAll(newIncompatiblePairs);

        // prune unneeded rotamers and rotamer pairs
        RotamerSpace.pruneIncompatibleRotamers(rotamerSpace, combinedIncompatiblePairs);

        // return new object
        //System.out.println("Finished pruning excess information.");
        return new DEEDoubles(rotamerSelfEnergies, newRotamerInteractionEnergies, rotamerSpace, combinedIncompatiblePairs);
    }

    /**
     * Makes a DEESingles object from this DEEDoubles object.  Basically just copies fields.
     */
    public DEESingles getDEESingles()
    {
        DEEDoubles newDoubles = this.prune();
        int numberOfRotamers = 0;
        for (List<Rotamer> list : rotamerSpace)
            numberOfRotamers += list.size();
        int incompatibleSize = newDoubles.oldIncompatiblePairs.size();
        //System.out.printf("New DEESingles object created with %d rotamers and %d incompatible pairs.\n", numberOfRotamers, incompatibleSize);
        return new DEESingles(rotamerSelfEnergies, newDoubles.rotamerInteractionEnergies, rotamerSpace, newDoubles.oldIncompatiblePairs);
    }

    /**
     * Calculates the epsilon value for two rotamers.  Epsilon(i_r, j_s) is defined as E(i_r) + E(j_s) + 
     * E(i_r,j_s), where i_r is rotamer r at position i and j_s is rotamer j at position j.
     * i should not equal j, but we do not check this.
     * @return the epsilon value, or null if an interesting rotamer or an incompatible pair is involed
     */
    public Double getEpsilon(Rotamer r, Rotamer s)
    {
        Double r_self = rotamerSelfEnergies.get(r);
        if ( r_self == null )
            return null;
        Double s_self = rotamerSelfEnergies.get(s);
        if ( s_self == null )
            return null;
        Double rs_interaction = rotamerInteractionEnergies.get(new RotamerPair(r,s));
        if ( rs_interaction == null )
            return null;
        return r_self + s_self + rs_interaction;
    }

    /** overload */
    public Double getEpsilon(RotamerPair pair)
    {
        Double r_self = rotamerSelfEnergies.get(pair.rotamer1);
        if ( r_self == null )
            return null;
        Double s_self = rotamerSelfEnergies.get(pair.rotamer2);
        if ( s_self == null )
            return null;
        Double rs_interaction = rotamerInteractionEnergies.get(pair);
        if ( rs_interaction == null )
            return null;
        return r_self + s_self + rs_interaction;
    }

    /**
     * Calculates the epsilon value for a rotamer pair (i_r, j_s) and another rotamer (k_t).  The energy is defined as
     * E(i_r,k_t) + E(j_s,k_t).  i, j, and k should not equal each other but we do not check this.
     * @return the epsilon value, or null if an interesting or incompatible pair is involved.
     */
    public Double getEpsilon(Rotamer r, Rotamer s, Rotamer t)
    {
        Double rt_interaction = rotamerInteractionEnergies.get(new RotamerPair(r,t));
        if ( rt_interaction == null )
            return null;
        Double st_interaction = rotamerInteractionEnergies.get(new RotamerPair(s,t));
        if ( st_interaction == null )
            return null;
        return rt_interaction + st_interaction;
    }

    /** overload */
    public Double getEpsilon(RotamerPair pair_rs, Rotamer t)
    {
        return getEpsilon(pair_rs.rotamer1, pair_rs.rotamer2, t);
    }

    /**
     * Precomputes epsilon_min and epsilon_max for all pairs (i_r, j_s).  Epsilon_min is defined as epsilon(r,s) +
     * sum over k [ max on t ( epsilon(r,s,t) ) ], where t refers to rotamer t at position k, k != i,j.  Epsilon_max
     * is defined as epsilon(r,s) + sum over k [ min on t ( epsilon(r,s,t) ) ], where t refers to rotamer t at
     * position k, k != i,j.
     * This does not actually do the precomputation, but simply sets up the jobs and combines the results.
     */
    public void preCompute()
    {
        if ( e_min_map.size() > 0 || e_max_map.size() > 0 )
            throw new IllegalArgumentException("maps already computed!");
        
        // reset debug counter that tells us how many pairs have been eliminated
        //debugCount.set(0);

        // enumerate all pairs and create jobs
        //System.out.println("Precomputing epsilon_max and epsilon_min for all valid rotamer pairs.");
        List<Future<Result>> futures = new LinkedList<>();
        int jobSize = 100;
        int pairCount = 0;
        int jobCount = 0;
        List<RotamerPair> thesePairs = new ArrayList<>(jobSize);
        for ( RotamerPair pair : rotamerInteractionEnergies.keySet() )
            {
                // if there is a transition state present then the one-center energy is not meaningful
                if ( pair.rotamer1.description.indexOf("transition_state") > -1 ||
                     pair.rotamer2.description.indexOf("transition_state") > -1    )
                    continue;
                pairCount++;
                thesePairs.add(pair);
                if ( pairCount >= jobSize )
                    {
                        jobCount++;
                        pairCount = 0;
                        PrecomputationJob job = new PrecomputationJob(thesePairs);
                        Future<Result> f = GeneralThreadService.submit(job);
                        futures.add(f);
                        //System.out.printf("Created %d precomputation jobs...\r", jobCount);
                        thesePairs = new ArrayList<>(jobSize);
                    }
            }

        // deal with the last batch
        if ( thesePairs.size() > 0 )
            {
                PrecomputationJob job = new PrecomputationJob(thesePairs);
                Future<Result> f = GeneralThreadService.submit(job);
                futures.add(f);
                //System.out.printf("Created %d precomputation jobs...\r", jobCount);
            }

        // wait for jobs to finish
        //System.out.println();
        GeneralThreadService.silentWaitForFutures(futures);
        //System.out.print("                                                    ");
        //System.out.println("\nPrecomputation complete.");
 
    }

    /**
     * Calculates epsilon_min and epsilon_max for a given list of pairs.  Saves time by doing the calculations together.
     */
    public class PrecomputationJob implements WorkUnit
    {
        public static final long serialVersionUID = 1L;

        public final List<RotamerPair> pairs;

        public PrecomputationJob(List<RotamerPair> pairs)
        {
            this.pairs = pairs;
        }

        public Result call()
        {
            int size = rotamerSpace.size();
            for (RotamerPair pair : pairs)
                {
                    // get fields
                    Rotamer r = pair.rotamer1;
                    int i = r.sequenceIndex;
                    Rotamer s = pair.rotamer2;
                    int j = s.sequenceIndex;

                    // initialize epsilon_max and epsilon_min
                    double maxSum = getEpsilon(pair);
                    double minSum = maxSum;

                    // add the maximum and minimum contribution from each position k
                    for (int k = 0; k < size; k++)
                        {
                            if ( k == i || k == j )
                                continue;
                            double thisMaxSum = 0.0;
                            double thisMinSum = 0.0;
                            boolean first = true;

                            // rotamer t is at position k
                            for (Rotamer t : rotamerSpace.get(k))
                                {
                                    Double energy = getEpsilon(r,s,t);
                                    if ( energy == null )
                                        continue;
                                    if ( first )
                                        {
                                            thisMaxSum = energy;
                                            thisMinSum = energy;
                                            first = false;
                                        }
                                    else
                                        {
                                            if ( energy > thisMaxSum )
                                                thisMaxSum = energy;
                                            if ( energy < thisMinSum )
                                                thisMinSum = energy;
                                        }
                                }
                            maxSum += thisMaxSum;
                            minSum += thisMinSum;
                        }

                    // update the results
                    if ( e_max_map.containsKey(pair) || e_min_map.containsKey(pair) )
                        throw new IllegalArgumentException("duplicate key");
                    e_max_map.put(pair, maxSum);
                    e_min_map.put(pair, minSum);
                }

            // because the results are updated concurrently, there is no need to return a real result
            return null;
        }
    }

    /**
     * Performs magic bullet doubles.
     */
    public void performMagicBulletDoubles()
    {
        // precompute the epsilon_max and epsilon_min maps
        preCompute();
        
        // for all pairs of positions, set up magic bullet jobs
        int size = rotamerSpace.size();
        List<MagicBulletJob> allJobs = new ArrayList<>();
        int latchSize = 0;
        for (int i=0; i < size-1; i++ )
            {
                // no pairs can be generated from this position so skip
                if ( rotamerSpace.get(i).size() <= 1 )
                    continue;
                for (int j=i+1; j < size; j++)
                    {
                        // no pairs can be generated from this position so skip
                        if ( rotamerSpace.get(j).size() <= 1 )
                            continue;
                        latchSize++;
                    }
            }
        CountDownLatch countDownLatch = new CountDownLatch(latchSize);
        List<Future<Result>> futures = new ArrayList<Future<Result>>();
        for (int i=0; i < size-1; i++ )
            {
                // no pairs can be generated from this position so skip
                if ( rotamerSpace.get(i).size() <= 1 )
                    continue;
                for (int j=i+1; j < size; j++)
                    {
                        // no pairs can be generated from this position so skip
                        if ( rotamerSpace.get(j).size() <= 1 )
                            continue;
                        MagicBulletJob job = new MagicBulletJob(i,j,futures,countDownLatch);
                        allJobs.add(job);
                    }
            }

        // submit jobs
        for (MagicBulletJob job : allJobs)
            {
                Future<Result> f = GeneralThreadService.submit(job);
                synchronized (futures) {
                    futures.add(f);
                }
            }

        // wait for jobs to complete
        //System.out.println("latch start");
        while (true)
            {
                try
                    {
                        countDownLatch.await();
                    }
                catch (InterruptedException e)
                    {
                        continue;
                    }
                break;
            }
        //System.out.println("latch done");
        while (true)
            {
                synchronized(futures)
                    {
                        int numberDone = 0;
                        for (Future<Result> f : futures)
                        if ( f.isDone() )
                            numberDone++;
                        int queueSize = GeneralThreadService.queueSize();
                        //System.out.printf("%d / %d    queue: %d           \r", numberDone, futures.size(), queueSize);
                        if ( numberDone == futures.size() && queueSize == 0 )
                            break;
                    }
                GeneralThreadService.wait(50);
            }
        //System.out.println();

        // make a report
        //System.out.printf("\nEliminated %d incompatible pairs.\n", newIncompatiblePairs.size());
    }

    //public static final AtomicInteger debugCount = new AtomicInteger();

    /**
     * Determines the magic bullet rotamer pair (u,v)_mb at positions i and j, and then examines all other rotamer pairs
     * (r,s) at i and j, such that (u,v) != (r,s), to see if they can be eliminated.
     */
    public class MagicBulletJob implements WorkUnit
    {
        /** the rotamer indices to perform magic bullet singles at */
        public final int i, j;
        
        /** the list where we put all the futures so we don't run off before the jobs are done */
        private final List<Future<Result>> futures;
        
        /** so the parent job knows to finish */
        private final CountDownLatch countDownLatch;

        public MagicBulletJob(int i, int j, List<Future<Result>> futures, CountDownLatch countDownLatch)
        {
            this.i = i;
            this.j = j;
            this.futures = futures;
            this.countDownLatch = countDownLatch;
        }
        
        public Result call()
        {
            // get fields
            List<Rotamer> listI = rotamerSpace.get(i);
            List<Rotamer> listJ = rotamerSpace.get(j);
            int size = rotamerSpace.size();

            // determine the magic bullet for this pair
            Rotamer u_mb = null;
            Rotamer v_mb = null;
            boolean first = true;
            Double epsilon_max_mb = null;
            for (Rotamer u : listI)
                {
                    for (Rotamer v : listJ)
                        {
                            RotamerPair uv = new RotamerPair(u,v);
                            Double epsilon_max = e_max_map.get(uv);
                            if ( epsilon_max == null )
                                continue;
                            if ( first || epsilon_max < epsilon_max_mb )
                                {
                                    u_mb = u;
                                    v_mb = v;
                                    epsilon_max_mb = epsilon_max;
                                    if ( first )
                                        first = false;
                                }
                        }
                }

            // if there is no valid magic bullet for this pair, return
            if ( u_mb == null || v_mb == null || epsilon_max_mb == null )
                return null;

            // get the magic bullet interaction energy
            RotamerPair uv_mb = new RotamerPair(u_mb, v_mb);
            Double energy_uv_mb = getEpsilon(uv_mb);

            // use the magic bullet to perform first-order doubles on all other pairs
            for (Rotamer r : listI)
                {
                    MagicBulletMicroJob job = new MagicBulletMicroJob(i, j, r, listJ, uv_mb, energy_uv_mb);
                    Future<Result> f = GeneralThreadService.submit(job);
                    synchronized (futures)
                        {
                            futures.add(f);
                        }
                }
            countDownLatch.countDown();
            //System.out.println(countDownLatch.getCount());
            return null;
        }
    }

    public class MagicBulletMicroJob implements WorkUnit
    {
        public static final long serialVersionUID = 1L;
        
        public final int i, j;
        public final Rotamer r;
        public final List<Rotamer> listJ;
        public final RotamerPair uv_mb;
        public final Double energy_uv_mb;

        public MagicBulletMicroJob(int i, int j, Rotamer r, List<Rotamer> listJ, RotamerPair uv_mb, Double energy_uv_mb)
        {
            this.i = i;
            this.j = j;
            this.r = r;
            this.listJ = listJ;
            this.uv_mb = uv_mb;
            this.energy_uv_mb = energy_uv_mb;
        }

        public Result call()
        {
            int size = rotamerSpace.size();
            labelS:
            for (Rotamer s : listJ)
                {
                    RotamerPair rs = new RotamerPair(r,s);
                    
                    // don't do the magic bullet
                    if ( rs.equals(uv_mb) )
                        continue;

                    // calculate the required sums and minimum operators
                    Double energy_rs = getEpsilon(rs);
                    if ( energy_rs == null )
                        continue; // could happen if rs is incompatible or at least one of r or s is interesting

                    double sum = energy_rs - energy_uv_mb;

                    for (int k=0; k < size; k++)
                        {
                            // skip i and j
                            if ( k==i || k==j ) // should only compute interactions with i--j and i--k
                                continue;

                            List<Rotamer> listK = rotamerSpace.get(k);
                            double thisMinimum = 0.0;
                            boolean first = true;
                            for (Rotamer t : listK)
                                {
                                    // get energies
                                    Double energy_rst = getEpsilon(rs, t);
                                    if ( energy_rst == null )
                                        continue; // happens if (r,t) or (s,t) is incompatible; irrelevant because t cannot coexist with (r,s)
                                    Double energy_uvt = getEpsilon(uv_mb, t);
                                    if ( energy_uvt == null )
                                        continue labelS;

                                    // perform the minimum operation
                                    double thisEnergyDifference = energy_rst - energy_uvt;
                                    if ( first || thisEnergyDifference < thisMinimum )
                                        {
                                            thisMinimum = thisEnergyDifference;
                                            if ( first )
                                                first = false;
                                        }
                                }
                            sum += thisMinimum;
                        }
                    
                    // check to see if this can be eliminated
                    if ( sum > 0.0 )
                        {
                            //debugCount.getAndIncrement();
                            //System.out.printf("%d,%d: magic bullet elimination of %s-%s\n", i, j,
                            //                  rs.rotamer1.protoAminoAcid.r.aminoAcid, rs.rotamer2.protoAminoAcid.r.aminoAcid);
                            newIncompatiblePairs.add(rs);
                        }
                }
            return null;
        }
    }
} // end of class DEEDoubles
