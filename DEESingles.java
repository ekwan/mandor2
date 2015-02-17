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
 * Performs zero-order singles elimination for the DEE algorithm.
 */
public class DEESingles
{
    /** the self-energy of the sidechain and the interaction energy with the backbone, corrected by reference energy */
    public final Map<Rotamer,Double> rotamerSelfEnergies;

    /** the interaction energies between the two rotamers */
    public final Map<RotamerPair,Double> rotamerInteractionEnergies;

    /** the space of all Rotamers allowed at the various positions on the amino acid.  */
    public final List<List<Rotamer>> rotamerSpace;

    /** the incompatible pairs */
    public final Set<RotamerPair> incompatiblePairs;

    /** constructor */
    public DEESingles(Map<Rotamer,Double> rotamerSelfEnergies,
                          Map<RotamerPair,Double> rotamerInteractionEnergies,
                          List<List<Rotamer>> rotamerSpace,
                          Set<RotamerPair> incompatiblePairs)
    {
        this.rotamerSelfEnergies = rotamerSelfEnergies;
        this.rotamerInteractionEnergies = rotamerInteractionEnergies;
        this.rotamerSpace = rotamerSpace;
        this.incompatiblePairs = incompatiblePairs;
    }

    public enum JobType
    {
        SPLIT_FIRST_ORDER, MB_SPLIT_FIRST_ORDER;
    }

    /**
     * Performs the elimination process.
     * @return the new rotamer space that has various rotamers eliminated
     */
    public List<List<Rotamer>> eliminate(JobType jobType)
    {
        // create jobs
        List<Future<Result>> futures = new LinkedList<>();
        for (int i=0; i < rotamerSpace.size(); i++)
            {
                // don't bother trying to eliminate anything at positions where there isn't anything to eliminate
                if ( rotamerSpace.get(i).size() <= 1 )
                    continue;

                WorkUnit job = null;
                for (int j = 0; j < rotamerSpace.get(i).size(); j++) {
                    Rotamer r = rotamerSpace.get(i).get(j);
                    if (jobType == JobType.SPLIT_FIRST_ORDER) job = new SplitFirstOrderSinglesJob(i,r);
                    else job = new MagicBulletSplitFirstOrderSinglesJob(i,r,j);
                    Future<Result> f = GeneralThreadService.submit(job);
                    futures.add(f);
                }
            }

        // wait for jobs to finish
        GeneralThreadService.silentWaitForFutures(futures);

        // combine results
        HashSet<Rotamer> allEliminated = new HashSet<>();
        for (Future<Result> f : futures)
            {
                SinglesResult result = null;
                try
                    {
                        result = (SinglesResult)f.get();
                    }
                catch (Exception e)
                    {
                        e.printStackTrace();
                    }
                allEliminated.addAll(result.eliminated);
            }

        int numberPruned = 0;
        List<List<Rotamer>> prunedRotamerSpace = new ArrayList<>();
        for (List<Rotamer> oldRotamerList : rotamerSpace)
            {
                ArrayList<Rotamer> newRotamerList = new ArrayList<>();
                for (Rotamer r : oldRotamerList)
                    {
                        if ( ! allEliminated.contains(r) )
                            newRotamerList.add(r);
                        else
                            numberPruned++;
                    }
                prunedRotamerSpace.add(newRotamerList);
            }
        //System.out.println(numberPruned + " rotamers have been pruned.");
        return prunedRotamerSpace;
    }

    /**
     * Produce a list of rotamers at a given position on the peptide that can be 
     * eliminated through split 1-order singles comparisons.  
     */
    public class SplitFirstOrderSinglesJob implements WorkUnit
    {
        /** for serialization */
        public static final long serialVersionUID = 1L;

        /** the position on the peptide where we are doing DEE */
        public final int i;

        /** the Rotamer we're trying to eliminate from position i.  */
        public final Rotamer r;

        /** 
         * Construct a split first-order singles job that works on position i of the rotamer space.
         * @param i the position of the peptide sequence where elimination is to occur
         */
        public SplitFirstOrderSinglesJob (int i, Rotamer r) {
            this.i = i;
            this.r = r;
        }

        /**
         * Does split 1-order singles DEE at position i, Rotamer r.  
         * @return a Result containing r and a HashSet containing r or nothing.  
         * if it contains r, then r is to be eliminated.  
         */
        public SinglesResult call()
        {
            boolean[] useless = new boolean[rotamerSpace.get(i).size()];
            double[][] energyMins = new double[rotamerSpace.size()][rotamerSpace.get(i).size()];
            // iterate through all other rotamers at this position
            MLoop:
            for (int m = 0; m < rotamerSpace.get(i).size(); m++) {
                Rotamer t = rotamerSpace.get(i).get(m);
                if (t.equals(r)) continue;
                // iterate through all other positions
                for (int j = 0; j < rotamerSpace.size(); j++) {
                    if (j==i) continue; 
                    double min = 0.0;
                    boolean first = true;
                    for (Rotamer s : rotamerSpace.get(j)) {
                        Double energy1 = rotamerInteractionEnergies.get(new RotamerPair(r,s));
                        Double energy2 = rotamerInteractionEnergies.get(new RotamerPair(t,s));
                        if (energy2==null) {
                            useless[m] = true;
                            continue MLoop;
                        }
                        if (energy1==null) continue;
                        double energy = energy1-energy2;
                        if (first) {
                            min = energy;
                            first = false;
                        } else {
                            if (energy < min) min = energy;
                        }
                    }
                    energyMins[j][m] = min;
                }
            }

            // we've got all the energy minima
            // now iterate through all k except i
            split:
            for (int k = 0; k < rotamerSpace.size(); k++) {
                if (k==i||rotamerSpace.get(k).size()==0) continue;
                // if we can get this to be true for all rotamers at k, then we eliminate r
                boolean[] elimk = new boolean[rotamerSpace.get(k).size()];
                comparand:
                for (int m = 0; m < rotamerSpace.get(i).size(); m++) {
                    if (useless[m]) continue;
                    Rotamer t = rotamerSpace.get(i).get(m);
                    if (t.equals(r)) continue;
                    double energy = rotamerSelfEnergies.get(r)-rotamerSelfEnergies.get(t);
                    for (int j = 0; j < rotamerSpace.size(); j++)
                        if ((j!=i)&&(j!=k)) energy += energyMins[j][m];
                    for (int n = 0; n < rotamerSpace.get(k).size(); n++) {
                        Rotamer kv = rotamerSpace.get(k).get(n);
                        Double energy1 = rotamerInteractionEnergies.get(new RotamerPair(r,kv));
                        Double energy2 = rotamerInteractionEnergies.get(new RotamerPair(t,kv));
                        if (energy1==null) {
                            elimk[n] = true;
                            continue;
                        }
                        if (energy2==null) continue comparand;
                        double extraEnergy = energy1-energy2;
                        if (energy + extraEnergy > 0.0) elimk[n] = true;
                    }
                }
                for (int n = 0; n < rotamerSpace.get(k).size(); n++) {
                    Rotamer kv = rotamerSpace.get(k).get(n);
                    // we found one partition in the split where we can't eliminate r
                    if (!elimk[n]) continue split;
                }
                // we can eliminate r at every spot in the partition
                HashSet<Rotamer> eliminated = new HashSet<Rotamer>();
                eliminated.add(r);
                //System.out.printf("Eliminated %s %s at position %d, splitting at position %d.                    \n", r.protoAminoAcid.r.aminoAcid, r.chis, i, k);
                return new SinglesResult(i,eliminated);
            }

            return new SinglesResult(i, new HashSet<Rotamer>());
        }
    } // end of class SplitFirstOrderSinglesJob

    /**
     * Produce a list of rotamers at a given position on the peptide that can be 
     * eliminated through split 1-order singles comparisons.  
     */
    public class MagicBulletSplitFirstOrderSinglesJob implements WorkUnit
    {
        /** for serialization */
        public static final long serialVersionUID = 1L;

        /** the position on the peptide where we are doing DEE */
        public final int i;

        /** the Rotamer we're trying to eliminate from position i.  */
        public final Rotamer r;

        /** the position of r in the list at peptide position i.  */
        public final int p;

        /** 
         * Construct a split first-order singles job that works on position i of the rotamer space.
         * @param i the position of the peptide sequence where elimination is to occur
         */
        public MagicBulletSplitFirstOrderSinglesJob (int i, Rotamer r, int p) {
            this.i = i;
            this.r = r;
            this.p = p;
        }

        /**
         * Does split 1-order singles DEE at position i, Rotamer r.  
         * @return a Result containing r and a HashSet containing r or nothing.  
         * if it contains r, then r is to be eliminated.  
         */
        public SinglesResult call()
        {
            double[][] energyMins = new double[rotamerSpace.size()][rotamerSpace.get(i).size()];
            boolean[] unusable = new boolean[rotamerSpace.get(i).size()];
            // iterate through all other rotamers at this position
            for (int m = 0; m < rotamerSpace.get(i).size(); m++) {
                if (m==p) {
                    unusable[m] = true;
                    continue;
                }
                Rotamer t = rotamerSpace.get(i).get(m);
                // iterate through all other positions
                for (int j = 0; j < rotamerSpace.size(); j++) {
                    if (j==i) continue;
                    double min = 0.0;
                    boolean first = true;
                    for (Rotamer s : rotamerSpace.get(j)) {
                        Double energy1 = rotamerInteractionEnergies.get(new RotamerPair(r,s));
                        Double energy2 = rotamerInteractionEnergies.get(new RotamerPair(t,s));
                        if (energy1==null) continue;
                        if (energy2==null) {
                            unusable[m] = true;
                            continue;
                        }
                        double energy = energy1-energy2;
                        if (first) {
                            min = energy;
                            first = false;
                        } else {
                            if (energy < min) min = energy;
                        }
                    }
                    energyMins[j][m] = min;
                }
            }

            // use the energyMins to find the best two splitting positions k1 and k2
            int k1 = -1;
            int k2 = -1;
            double min1 = 0.0;
            double min2 = 0.0;
            int count = 0;
            for (int k = 0; k < rotamerSpace.size(); k++) {
                if (k==i||rotamerSpace.get(k).size()==0) continue;
                count++;
                double min = 0.0;
                // iterate through all competitors t at i
                boolean first = true;
                for (int m = 0; m < energyMins[0].length; m++) {
                    if (m==p) continue;
                    if (first) {
                        min = energyMins[k][m];
                        first = false;
                    } else {
                        if (energyMins[k][m] < min) min = energyMins[k][m];
                    }
                }
                if (count==1) {
                    min1 = min;
                    k1 = k;
                } else if (count == 2) {
                    if (min < min1) {
                        // bump the current best value down
                        min2 = min1;
                        min1 = min;
                        k2 = k1;
                        k1 = k;
                    } else {
                        // put the min we just found in the second value
                        min2 = min;
                        k2 = k;
                    }
                } else {
                    if (min < min2) {
                        if (min < min1) {
                            // bump the current best value down
                            min2 = min1;
                            min1 = min;
                            k2 = k1;
                            k1 = k;
                        } else {
                            // put the min we just found in the second value
                            min2 = min;
                            k2 = k;
                        }
                    }
                }
            }

            // we've got all the energy minima
            // we've got k1 and k2
            // check every element at k1 and k2 to see if r can be eliminated there
            for (int x = 0; x < rotamerSpace.get(k1).size(); x++) {
                for (int y = 0; y < rotamerSpace.get(k2).size(); y++) {
                    // flag tells us if r can be eliminated at rotamer x of position k1 and rotamer y of position k2
                    boolean elim = false;
                    // now test all competitor rotamers at i
                    for (int m = 0; m < rotamerSpace.get(i).size(); m++) {
                        if (unusable[m]) continue;
                        Rotamer t = rotamerSpace.get(i).get(m);
                        elim = false;
                        // start with self-energy difference between r and t
                        double energy = rotamerSelfEnergies.get(r)-rotamerSelfEnergies.get(t);
                        // add intereaction energy difference with kv, which is at k1
                        Rotamer kv = rotamerSpace.get(k1).get(x);
                        Double energy11 = rotamerInteractionEnergies.get(new RotamerPair(r,kv));
                        Double energy12 = rotamerInteractionEnergies.get(new RotamerPair(t,kv));
                        if (energy11==null) {
                            elim = true;
                            break;
                        }
                        energy += energy11-energy12;
                        // add intereaction energy difference with ku, which is at k2
                        Rotamer ku = rotamerSpace.get(k2).get(y);
                        Double energy21 = rotamerInteractionEnergies.get(new RotamerPair(r,ku));
                        Double energy22 = rotamerInteractionEnergies.get(new RotamerPair(t,ku));
                        if (energy21==null) {
                            elim = true;
                            break;
                        }
                        energy += energy21-energy22;
                        for (int j = 0; j < rotamerSpace.size(); j++)
                            if ((j!=i)&&(j!=k1)&&(j!=k2)) energy += energyMins[j][m];
                        if (energy > 0.0) {
                            elim = true;
                            break;
                        }
                    }
                    if (!elim) return new SinglesResult(i, new HashSet<Rotamer>());
                }
            }
            // we can eliminate r at every spot in the partition
            HashSet<Rotamer> eliminated = new HashSet<Rotamer>();
            eliminated.add(r);
            //System.out.printf("Eliminated %s %s at position %d, splitting at magic bullet positions %d and %d.                    \n", r.protoAminoAcid.r.aminoAcid, r.chis, i, k1, k2);
            return new SinglesResult(i,eliminated);
        }
    } // end of class MagicBulletSplitFirstOrderSinglesJob

    /**
     * The result of a SinglesJob.
     */
    public class SinglesResult implements Result
    {
        /** for serialization */
        public static final long serialVersionUID = 1L;

        /** the position on the peptide where we are doing DEE */
        public final int sequenceIndex;

        /** the rotamers that should be eliminated from the rotamer space*/
        public final HashSet<Rotamer> eliminated;

        public SinglesResult(int sequenceIndex, HashSet<Rotamer> eliminated) {
            this.sequenceIndex = sequenceIndex;
            this.eliminated = eliminated;
        }
    } // end of class SinglesResult
} // end of class DEESingles
