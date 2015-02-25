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
 * This class computes the one- and two-center energy terms for a set of rotamers.
 * The backbone is defined as the set of atoms common to all conformations.  
 * Because this is designed to be used with FixedSequenceRotamerSpace, 
 * we will never go from proline to non-proline, or vice-versa.  Therefore, the
 * backbone HNs are always present and will not be placed in the sidechains.
 *
 * Reference energies are unnecessary here because we are dealing with a fixed composition.
 */
public class AStarEnergyCalculator
{
    /** the self-energy of the sidechain and the interaction energy with the backbone */
    public final Map<Rotamer,Double> rotamerSelfEnergies;

    /** the interaction energies between the two rotamers */
    public final Map<RotamerPair,Double> rotamerInteractionEnergies;

    /** the backbone self energy */
    public final double backboneEnergy;

    /** constructor */
    private AStarEnergyCalculator(Map<Rotamer,Double> rotamerSelfEnergies, Map<RotamerPair,Double> rotamerInteractionEnergies, double backboneEnergy)
    {
        this.rotamerSelfEnergies = rotamerSelfEnergies;
        this.rotamerInteractionEnergies = rotamerInteractionEnergies;
        this.backboneEnergy = backboneEnergy;
    }

    /** parallelization is turned on by default */
    public static AStarEnergyCalculator analyze(FixedSequenceRotamerSpace fixedSequenceRotamerSpace)
    {
        return analyze(fixedSequenceRotamerSpace, true);
    }

    /**
     * Factory method to create a AStarEnergyCalculator.  Calculates the single and pairwise energies of all the rotamers
     * in the RotamerPacker.  Energies on a truncated OPLS forcefield containing only charge and vdw interactions.
     * @param fixedSequenceRotamerSpace the rotamer space to analyze
     * @param parallelize whether to parallelize the calculation
     * @return the calculator
     */
    public static AStarEnergyCalculator analyze(FixedSequenceRotamerSpace fixedSequenceRotamerSpace, boolean parallelize)
    {
        // get some information
        Peptide peptide = fixedSequenceRotamerSpace.peptide;
        List<List<Rotamer>> rotamerSpace = fixedSequenceRotamerSpace.rotamerSpace;
        Set<RotamerPair> incompatiblePairs = fixedSequenceRotamerSpace.incompatiblePairs;

        // estimate the number of rotamers and rotamer pairs
        int totalRotamers = RotamerSpace.countRotamers(rotamerSpace);;
        int estimated = ( ( totalRotamers * (totalRotamers - 1) ) / 2 ) - incompatiblePairs.size();

        // populate fields
        ConcurrentHashMap<Rotamer,Double> selfEnergiesMap = new ConcurrentHashMap<Rotamer,Double>(totalRotamers, 0.75F, 16);
        ConcurrentHashMap<RotamerPair,Double> interactionEnergiesMap = new ConcurrentHashMap<RotamerPair,Double>(estimated, 0.75F, 16);
        AtomicDouble backboneEnergy = new AtomicDouble();
        
        // create energy jobs
        List<Future<Result>> futures = createEnergyJobs(rotamerSpace, peptide, incompatiblePairs, selfEnergiesMap, interactionEnergiesMap, backboneEnergy, parallelize);
        GeneralThreadService.silentWaitForFutures(futures);

        // return the object
        return new AStarEnergyCalculator(selfEnergiesMap, interactionEnergiesMap, backboneEnergy.get());
    }

    /**
     * Debugging method that predicts the energy of a peptide constructed from a list of rotamers.
     * @param backboneEnergy the energy of the backbone
     * @param rotamerSelfEnergies the one-center energies
     * @param rotamerInteractionEnergies the two-center energies
     * @param rotamerList the rotamers
     * @return the predicted energy in kcal
     */
    public static double predictEnergy(double backboneEnergy, Map<Rotamer,Double> rotamerSelfEnergies, Map<RotamerPair,Double> rotamerInteractionEnergies, List<Rotamer> rotamerList)
    {
        // add the backbone energy
        double predictedEnergy = backboneEnergy;
        //System.out.printf("backbone : %.4f\n", backboneEnergy);

        // add the energies of all the rotamers
        for (Rotamer r : rotamerList)
            {
                if ( ! rotamerSelfEnergies.containsKey(r) )
                    throw new IllegalArgumentException("missing energy for rotamer " + r.toString() + "  index = " + rotamerList.indexOf(r));
                predictedEnergy += rotamerSelfEnergies.get(r);
                //System.out.printf("%s (%d) : %.4f\n", r.protoAminoAcid.r.description, r.sequenceIndex, rotamerSelfEnergies.get(r));
            }

        // add the interaction energies
        for (int i=0; i < rotamerList.size(); i++)
            {
                Rotamer rotamer1 = rotamerList.get(i);
                for (int j=i+1; j < rotamerList.size(); j++)
                    {
                        Rotamer rotamer2 = rotamerList.get(j);
                        RotamerPair pair = new RotamerPair(rotamer1, rotamer2);
                        if ( ! rotamerInteractionEnergies.containsKey(pair) )
                            throw new IllegalArgumentException("missing energy for rotamer pair: " + pair.toString());
                        Double energy = rotamerInteractionEnergies.get(pair);
                        predictedEnergy += energy;
                        //System.out.printf("%s (%d) / rotamer %s (%d) : %.4f\n", rotamer1.protoAminoAcid.r.description,
                        //                  rotamer1.sequenceIndex, rotamer2.protoAminoAcid.r.description, rotamer2.sequenceIndex, energy);
                    }
            }
        return predictedEnergy;
    }

    /**
     * Creates a matrix of rotamer energies for a particular peptide.  Diagonal entries are rotamer self-energies,
     * which are defined as the interactions inside the rotamer as well as rotamer-backbone interactions.  The
     * entries are arranged in sequence order (0, 1, ..., n-1).  For example, entry (0,1) (or (1,0), since the matrix
     * is symmetric) corresponds to the interaction energy between the rotamer at position 0 and the rotamer at position 1.
     * The backbone is treated as the n-th residue, such that the (n,n) entry is the backbone self-energy.
     *
     * This variant does not put backbone HNs in the rotamer sidechains.  Rotamers will be placed at all non-hairpin positions
     * so every position must have at least one entry in its rotamer space.
     *
     * @param peptide the peptide that contains the rotamers
     * @param interactions the interactions in this peptide
     * @return the matrix of rotamer/backbone - rotamer/backbone interaction energies
     */
    public static Double[][] getRotamerEnergyMatrix(Peptide peptide, List<Interaction> interactions)
    {
        // initialize matrix where the results will go
        // the backbone is treated as the n+1-th residue,
        // where n is the number of residues
        int numberOfResidues = peptide.sequence.size();
        double[][] energyMatrix = new double[numberOfResidues+1][numberOfResidues+1];

        // get all rotamer atoms
        Set<Atom> allRotamerAtoms = new HashSet<>();
        List<Set<Atom>> rotamerAtoms = new ArrayList<>(numberOfResidues);
        for (int i=0; i < peptide.sequence.size(); i++)
            {
                Residue residue = peptide.sequence.get(i);
                Set<Atom> set = new HashSet<>();
                if ( !residue.isHairpin )
                    set = RotamerFactory.getSidechainAtoms(peptide, residue, false);
                rotamerAtoms.add(set);
                allRotamerAtoms.addAll(set);
            }
        
        // get backbone atoms
        Set<Atom> backboneAtoms = new HashSet<>(peptide.contents);
        backboneAtoms.removeAll(allRotamerAtoms);

        // classify the interactions
        for (Interaction interaction : interactions)
            {
                List<Atom> backboneAtomsList = ImmutableList.copyOf(interaction.atoms);

                // figure out which groups this interaction belongs to
                boolean inBackbone = false;
                boolean inRotamer1 = false;
                int rotamer1index  = -1;
                boolean inRotamer2 = false;
                int rotamer2index  = -1;

                for (Atom a : interaction.atoms)
                    {
                        if ( backboneAtoms.contains(a) )
                            inBackbone = true;
                        else
                            {
                                for (int rotamerIndex = 0; rotamerIndex < rotamerAtoms.size(); rotamerIndex++)
                                    {
                                        Set<Atom> thisRotamerAtoms = rotamerAtoms.get(rotamerIndex);
                                        if ( thisRotamerAtoms.contains(a) )
                                            {
                                                if ( rotamerIndex == rotamer1index || rotamerIndex == rotamer2index )
                                                    break;
                                                else if ( rotamer1index == -1 )
                                                    {
                                                        inRotamer1 = true;
                                                        rotamer1index = rotamerIndex;
                                                        break;
                                                    }
                                                else if ( rotamer2index == -1 )
                                                    {
                                                        inRotamer2 = true;
                                                        rotamer2index = rotamerIndex;
                                                        break;
                                                    }
                                                else
                                                    throw new IllegalArgumentException("interaction cannot be in three rotamers");
                                            }
                                    }
                            }
                    }

                // put the energy in the correct bin
                double interactionEnergy = interaction.interactionEnergy;
                if      (  inBackbone && !inRotamer1 && !inRotamer2 )
                    {
                        //System.out.println("backbone only");
                        energyMatrix[numberOfResidues][numberOfResidues] += interactionEnergy;
                        /*
                        backboneInteractions++;
                        backboneEnergy += interactionEnergy;
                        
                        int number1 = backboneList.indexOf(backboneAtomsList.get(0))+1;
                        int number2 = backboneList.indexOf(backboneAtomsList.get(1))+1;
                        if ( number1 < number2 )
                            backboneDescriptionList.add(String.format("\n%d-%d  :  %s", number1, number2, interaction.description));
                        else
                            backboneDescriptionList.add(String.format("\n%d-%d  :  %s", number2, number1, interaction.description));
                        */
                    }
                else if (  inBackbone &&  inRotamer1 && !inRotamer2 )
                    {
                        //System.out.println("position " + rotamer1index + " backbone");
                        energyMatrix[rotamer1index][numberOfResidues] += interactionEnergy;
                        energyMatrix[numberOfResidues][rotamer1index] += interactionEnergy;
                    }
                else if (  inBackbone && !inRotamer1 &&  inRotamer2 )
                    {
                        //System.out.println("position " + rotamer2index + " backbone");
                        energyMatrix[rotamer2index][numberOfResidues] += interactionEnergy;
                        energyMatrix[numberOfResidues][rotamer2index] += interactionEnergy;
                    }
                else if ( !inBackbone &&  inRotamer1 && !inRotamer2 )
                    {
                        //System.out.println("position " + rotamer1index + " only");
                        energyMatrix[rotamer1index][rotamer1index] += interactionEnergy;
                    }
                else if ( !inBackbone && !inRotamer1 &&  inRotamer2 )
                    {
                        //System.out.println("position " + rotamer2index + " only");
                        energyMatrix[rotamer2index][rotamer2index] += interactionEnergy;
                    }
                else if ( !inBackbone &&  inRotamer1 &&  inRotamer2 )
                    {
                        //System.out.println("position " + rotamer1index + ", " + rotamer2index + " interaction");
                        energyMatrix[rotamer1index][rotamer2index] += interactionEnergy;
                        energyMatrix[rotamer2index][rotamer1index] += interactionEnergy;
                    }
                else
                    throw new IllegalArgumentException("error assigning interaction");
            }
        
        /*System.out.printf("[%3d] %d backbone self-interactions, %d backbone atoms, %.4f kcal\n", count, backboneInteractions, backboneAtoms.size(), backboneEnergy);

        String debugFilename3 = String.format("test_peptides/backbone_%05d.txt", count);
        
        String backboneDescription = "";
        Collections.sort(backboneDescriptionList);
        for (String s : backboneDescriptionList)
            backboneDescription += s;
        InputFileFormat.writeStringToDisk(backboneDescription, debugFilename3);
        */

        // return result
        Double[][] resultMatrix = new Double[numberOfResidues+1][numberOfResidues+1];
        for (int i=0; i < numberOfResidues+1; i++)
            {
                for (int j=0; j < numberOfResidues+1; j++)
                    {
                        resultMatrix[i][j] = Double.valueOf(energyMatrix[i][j]);
                    }
            }

        return resultMatrix;
    }

    /**
     * Creates all the RotamerEnergyJobs we need to get rotamer and rotamer pair energies
     * for an entire rotamer space.
     * @param rotamerSpace all the rotamers to analyze
     * @param startingPeptide contains the peptide backbone
     * @param incompatiblePairs all the incompatible pairs of rotamers we don't need interaction energies for
     * @param selfEnergiesMap the map to concurrently update for rotamer single energies
     * @param interactionEnergiesMap the map to concurrently update for rotamer pair energies
     * @param backboneEnergy the backbone energy
     * @param parallelize whether to do things in parallel
     * @return the results of all the calculations (dummy objects in this case)
     */
    public static List<Future<Result>> createEnergyJobs(List<List<Rotamer>> rotamerSpace, Peptide startingPeptide, Set<RotamerPair> incompatiblePairs,
                                                        ConcurrentHashMap<Rotamer,Double> selfEnergiesMap, ConcurrentHashMap<RotamerPair,Double> interactionEnergiesMap,
                                                        AtomicDouble backboneEnergy, boolean parallelize)
    {
        // create rotamer jobs
        // we minimize the amount of work by going across each row of rotamers
        // that is, we take the first rotamer at each position and make a peptide
        // then we go down the list, skipping positions where there aren't enough rotamers
        // incompatibles should not be an issue because those interactions will be ignored
        List<Future<Result>> futures = new ArrayList<>(); 
        int maxIndex = 0;
        for (List<Rotamer> list : rotamerSpace)
            {
                if ( list.size() - 1 > maxIndex )
                    maxIndex = list.size() -1;
            }

        for (int i=0; i <= maxIndex; i++)
            {
                List<Rotamer> rotamers = new ArrayList<>();
                for (List<Rotamer> list : rotamerSpace)    
                    {
                        if ( i < list.size() )
                            rotamers.add( list.get(i) );
                    }
                RotamerEnergyJob job = new RotamerEnergyJob(startingPeptide, rotamers, selfEnergiesMap, backboneEnergy);
                if ( parallelize )
                    {
                        //System.out.println("submit 2");
                        Future<Result> f = GeneralThreadService.submit(job);
                        futures.add(f);
                    }
                else
                    job.call();
            }

        // create rotamer pair jobs
        // we save time by just computing the interactions between the rotamers without the backbone present
        PairIterator iterator = new PairIterator(rotamerSpace, 100);
        while (iterator.hasNext())
            {
                LinkedListMultimap<Rotamer,Rotamer> thisBatch = iterator.next();
                RotamerPairEnergyJob job = new RotamerPairEnergyJob(startingPeptide, incompatiblePairs, thisBatch, interactionEnergiesMap);
                if ( parallelize )
                    {
                        //System.out.println("submit 3");
                        Future<Result> f = GeneralThreadService.submit(job);
                        futures.add(f);
                    }
                else
                    job.call();
            }
        return futures;
    }

    /**
     * Given a rotamer space, create all the possible pairs of rotamers without
     * any duplicates or intra-residue pairs.
     */
    public static class PairIterator
    {
        /** the rotamer space to traverse */
        private final List<List<Rotamer>> rotamerSpace;

        /** whether the iterator has anything to give next */
        private boolean hasNext;

        /** how many rotamer pairs to return at a time */
        private final int batchSize;

        /**
         * the current pair we are producing
         */
        private IndexPair currentPosition;

        /** constructor */
        public PairIterator(List<List<Rotamer>> rotamerSpace, int batchSize)
        {
            this.rotamerSpace = rotamerSpace;
            this.batchSize = batchSize;

            // find the first residue that has a non-zero number of rotamers
            int index1 = -1;
            int index3 = -1;
            for (List<Rotamer> list : rotamerSpace)
                {
                    if ( index1 == -1 && list.size() > 0 )
                        index1 = rotamerSpace.indexOf(list);
                    else if ( index3 == -1 && list.size() > 0 )
                        index3 = rotamerSpace.indexOf(list);
                }
            if ( index1 == -1 || index3 == -1 )
                throw new IllegalArgumentException("nothing to start the iterator with");
            hasNext = true;
            currentPosition = new IndexPair(index1,0,index3,0);
        }

        /**
         * Returns the next list we can use as the source for a rotamer pair.
         * @param currentFromList the current source
         * @return the next list
         */
        private List<Rotamer> nextFromList(List<Rotamer> currentFromList)
        {
            List<Rotamer> nextFromList = null;
            for (int i=rotamerSpace.indexOf(currentFromList)+1; i < rotamerSpace.size() - 1; i++)
                {
                    List<Rotamer> list = rotamerSpace.get(i);
                    if ( list.size() > 0 )
                        {
                            nextFromList = list;
                            break;
                        }
                }
            return nextFromList;
        }

        /**
         * Returns the next list we can use as the target for a rotamer pair.
         * @param currentToList the current target
         * @return the next list
         */
        private List<Rotamer> nextToList(List<Rotamer> currentToList)
        {
            List<Rotamer> nextToList = null;
            for (int i=rotamerSpace.indexOf(currentToList)+1; i < rotamerSpace.size(); i++)
                {
                    List<Rotamer> list = rotamerSpace.get(i);
                    if ( list.size() > 0 )
                        {
                            nextToList = list;
                            break;
                        }
                }
            return nextToList;
        }

        /**
         * Whether we are finished with this iterator.
         * @return false if we are done
         */
        public boolean hasNext()
        {
            return hasNext;
        }

        /**
         * Return the next batch of rotamer pairs.  The currentPosition will be set at
         * the first thing of the next map to return on the next call of next().
         * Access the resulting pairs by iterating over the entries() view of the multimap.
         * @return the next batch of rotamer pairs
         */
        public LinkedListMultimap<Rotamer,Rotamer> next()
        {
            // if we are finished, throw an exception
            if ( !hasNext )
                throw new IllegalArgumentException("iterator is finished");
            // otherwise, make the next batch
            
            // the current position
            int index1 = currentPosition.outerIndex1;
            int index2 = currentPosition.innerIndex1;
            int index3 = currentPosition.outerIndex2;
            int index4 = currentPosition.innerIndex2;

            // create the next batch
            List<Rotamer> fromList = rotamerSpace.get(index1);
            Rotamer fromRotamer = fromList.get(index2);
            List<Rotamer> toList = rotamerSpace.get(index3);
            Rotamer toRotamer = toList.get(index4);

            // iterate through the current batch
            LinkedListMultimap<Rotamer,Rotamer> nextBatch = LinkedListMultimap.create();
            while ( nextBatch.size() < batchSize )
                {
                    // add the current pair
                    nextBatch.put(fromRotamer,toRotamer);
                    
                    // increment the current position
                    index4++;

                    // if we have run out for the current to list
                    if ( index4 > toList.size() - 1 )
                        {
                            // try to move onto the next to list
                            index4 = 0;
                            toList = nextToList(toList);

                            // if we have run out of to lists
                            if ( toList == null )
                                {
                                    // try to move onto the next entry in the from list
                                    index2++;

                                    // if we have run out for the current from list
                                    if ( index2 > fromList.size() - 1 )
                                        {
                                            // we have run out of rotamers in this from list, so
                                            // attempt to move to the next from list
                                            index2 = 0;
                                            fromList = nextFromList(fromList);

                                            if ( fromList != null )
                                                {
                                                    // get the first rotamer
                                                    fromRotamer = fromList.get(index2);
                                                    toList = nextToList(fromList);
                                                    if ( toList == null )
                                                        {
                                                            hasNext = false;
                                                            return nextBatch;
                                                        }
                                                }
                                            else
                                                {
                                                    // we have run out of from lists, so we are finished iterating
                                                    hasNext = false;
                                                    return nextBatch;
                                                }
                                        }
                                    else
                                        {
                                            // entries are still available in the from list, so
                                            // get the new from rotamer
                                            fromRotamer = fromList.get(index2);
                                            toList = nextToList(fromList);
                                            if ( toList == null )
                                                {
                                                    hasNext = false;
                                                    return nextBatch;
                                                }
                                        }
                                }
                        }

                    // set the next pair
                    toRotamer = toList.get(index4);
                }

            // set the current position
            index1 = rotamerSpace.indexOf(fromList);
            index3 = rotamerSpace.indexOf(toList);
            currentPosition = new IndexPair(index1,index2,index3,index4);

            // return the current next
            return nextBatch;
        }
    }

    /**
     * Calculates the self-energies of the specified peptides.
     */
    public static class RotamerEnergyJob implements WorkUnit
    {
        /** for serialization */
        public static final long serialVersionUID = 1L;

        /** the starting peptide */
        public final Peptide startingPeptide;

        /** the rotamers to obtain the energies of */
        public final List<Rotamer> rotamers;

        /** where to put the rotamer energies */
        public final ConcurrentHashMap<Rotamer,Double> map;

        /** where to put the backbone energy */
        public final AtomicDouble backboneEnergy;

        /**
         * Constructor.  We don't check that everything is consistent.  In particular, we don't check that rotamerToReconstitute
         * is compatible with rotamers and rotamerPairs.
         * @param startingPeptide the starting peptide, which might not have the correct rotamers
         * @param rotamers the rotamers to compute the energies of
         * @param map the map to concurrently update with single rotamer energies
         * @param backboneEnergy the backbone energy
         */
        public RotamerEnergyJob(Peptide startingPeptide, List<Rotamer> rotamers, ConcurrentHashMap<Rotamer,Double> map, AtomicDouble backboneEnergy)
        {
            this.startingPeptide = startingPeptide;
            this.rotamers = rotamers;
            this.map = map;
            this.backboneEnergy = backboneEnergy;
        }

        /**
         * Creates the necessary peptide, calculates all the interactions with tinker, classifies all interactions,
         * and then populates the requested lists of 
         */
        public Result call()
        {
            // create peptide
            Peptide peptide = Rotamer.reconstitute(startingPeptide, rotamers);

            // analyze interactions
            List<Interaction> interactions = OPLScalculator.getInteractions(peptide);

            // initialize maps that will contain the answers
            Map<Rotamer,Double> rotamerMap = new LinkedHashMap<>();

            // analyze all interactions in this peptide
            // note that we call a special getRotamerEnergyMatrix() method specific to this class,
            // which treats all backbone HNs as part of the backbone and not the rotamer, as is done
            // in DEEenergyCalculator
            Double[][] energyMatrix = getRotamerEnergyMatrix(peptide, interactions);

            // get rotamer energies
            for (Rotamer rotamer : rotamers)
                {
                    // get the rotamer energy without the reference energy
                    Double singleEnergy = Interaction.getRotamerEnergy(rotamer, energyMatrix, false); 
                    rotamerMap.put(rotamer,singleEnergy);
                }

            // get backbone energy
            double thisBackboneEnergy = Interaction.getBackboneEnergy(energyMatrix); 

            // concurrently update results
            map.putAll(rotamerMap);
            backboneEnergy.compareAndSet(0.0, thisBackboneEnergy);
            return null;
        }

        @Override
        public String toString()
        {
            return String.format("RotamerEnergyJob for %d rotamers", rotamers.size());
        }
    }

    /**
     * Calculates the interaction energies for the specified rotamers.
     */
    public static class RotamerPairEnergyJob implements WorkUnit
    {
        /** for serialization */
        public static final long serialVersionUID = 1L;

        public final Peptide peptide;

        /** to access the incompatible list */
        public final Set<RotamerPair> incompatiblePairs;

        /** the rotamer pairs to analyze */
        public final LinkedListMultimap<Rotamer,Rotamer> work;

        /** where to put the results (map from pairs of rotamers to their interaction energies) */
        public final ConcurrentHashMap<RotamerPair,Double> map;

        /**
         * Create a batch of rotamer pairs to compute the interaction energy of.
         */
        public RotamerPairEnergyJob(Peptide peptide, Set<RotamerPair> incompatiblePairs, LinkedListMultimap<Rotamer,Rotamer> work, ConcurrentHashMap<RotamerPair,Double> map)
        {
            this.peptide = peptide;
            this.incompatiblePairs = incompatiblePairs;
            this.work = work;
            this.map = map;
        }
        
        /**
         * Computes the rotamer interaction energies.
         */
        public Result call()
        {
            HashMap<RotamerPair,Double> tempMap = new HashMap<>();
            for ( Map.Entry<Rotamer,Rotamer> entry : work.entries() )
                {
                    Rotamer rotamer1 = entry.getKey();
                    Rotamer rotamer2 = entry.getValue();
                    RotamerPair pair = new RotamerPair(rotamer1, rotamer2);
                    if ( incompatiblePairs.contains(pair) )
                        continue;
                    
                    // calculate energies but don't include backbone HNs
                    double energy = OPLScalculator.getInteractionEnergy(rotamer1, null, rotamer2, null);
                    tempMap.put(pair,energy);
                }
            map.putAll(tempMap);
            return null;
        }
    }

    /** for debugging */
    public static void compareBackboneAtoms(Peptide p1, Peptide p2)
    {
        // get all rotamer atoms
        int numberOfResidues = p1.sequence.size();
        List<Set<Atom>> rotamerAtoms = new ArrayList<>(numberOfResidues);
        Set<Atom> allRotamerAtoms = new HashSet<>();
        for (Residue r : p1.sequence)
            {
                Set<Atom> atoms = new HashSet<>();
                if ( !r.isHairpin )
                    atoms = RotamerFactory.getSidechainAtoms(p1,r,false);
                rotamerAtoms.add(atoms);
                allRotamerAtoms.addAll(atoms);
            }

        // get backbone atoms
        Set<Atom> backboneAtoms = new HashSet<>();
        for (Atom a : p1.contents)
            {
                if ( allRotamerAtoms.contains(a) )
                    continue;
                backboneAtoms.add(a);
            }

        // get all rotamer atoms
        List<Set<Atom>> rotamerAtoms2 = new ArrayList<>(numberOfResidues);
        Set<Atom> allRotamerAtoms2 = new HashSet<>();
        for (Residue r : p2.sequence)
            {
                Set<Atom> atoms = new HashSet<>();
                if ( !r.isHairpin )
                    atoms = RotamerFactory.getSidechainAtoms(p2,r,false);
                rotamerAtoms2.add(atoms);
                allRotamerAtoms2.addAll(atoms);
            }

        // get backbone atoms
        Set<Atom> backboneAtoms2 = new HashSet<>();
        for (Atom a : p2.contents)
            {
                if ( allRotamerAtoms2.contains(a) )
                    continue;
                backboneAtoms2.add(a);
            }

        List<Atom> list1 = new ArrayList<>(backboneAtoms);
        List<Atom> list2 = new ArrayList<>(backboneAtoms2);
        Collections.sort(list1);
        Collections.sort(list2);
        System.out.println(list1.size());
        System.out.println(list2.size());
        for (int i=0; i < list1.size(); i++)
            {
                Atom a1 = list1.get(i);
                Atom a2 = list2.get(i);
                System.out.printf("%s : %s : %3d : %.4f\n", a1.toString(), a2.toString(), a1.type2-a2.type2, Molecule.getDistance(a1,a2)); 
            }
    }

    /** for testing */
    public static void main(String[] args)
    {
        // create a peptide
        DatabaseLoader.go();
        List<Peptide> sheets = BetaSheetGenerator.generateSheets(6, 5, 10000, 0.01);
        Peptide peptide = sheets.get(0);
        
        int sequenceLength = peptide.sequence.size();
        List<String> stringSequence = ImmutableList.of("glycine", "valine", "asparagine", "aspartate", "glutamine", "valine",
                                                       "histidine_hd", "isoleucine", "phenylalanine", "serine", "threonine", "standard_alanine");
        List<ProtoAminoAcid> protoAminoAcids = ProtoAminoAcidDatabase.getSpecificSequence(stringSequence);
        int tempJ = 0;
        for (int i=0; i < sequenceLength; i++)
            {
                Residue residue = peptide.sequence.get(i);
                if ( residue.isHairpin )
                    continue;
                ProtoAminoAcid protoAminoAcid = protoAminoAcids.get(tempJ);
                peptide = SidechainMutator.mutateSidechain(peptide, residue, protoAminoAcid);
                tempJ++; 
            }

        // create the A* iterator
        // note that the includeHN should be set to false if we want to do A* here
        new GaussianInputFile(peptide).write("test_peptides/original.gjf");
        FixedSequenceRotamerSpace fixedSequenceRotamerSpace = new FixedSequenceRotamerSpace(peptide, false);
        List<List<Rotamer>> rotamerSpace = fixedSequenceRotamerSpace.rotamerSpace;
        AStarEnergyCalculator calculator = AStarEnergyCalculator.analyze(fixedSequenceRotamerSpace);
        double backboneEnergy = calculator.backboneEnergy;
        Map<Rotamer,Double> rotamerSelfEnergies = calculator.rotamerSelfEnergies;
        Map<RotamerPair,Double> rotamerInteractionEnergies = calculator.rotamerInteractionEnergies; 
        RotamerIterator iterator = new RotamerIterator(rotamerSpace, rotamerSelfEnergies, rotamerInteractionEnergies, 1000);
        
        // perform A* iteration, checking to see if the predicted and actual energies are the same
        List<RotamerIterator.Node> solutions = iterator.iterate();
        int count = 0;
        for (RotamerIterator.Node node : solutions)
            {
                List<Rotamer> rotamers = node.rotamers;
                
                // get predicted energy
                double predictedEnergy = AStarEnergyCalculator.predictEnergy(backboneEnergy, rotamerSelfEnergies, rotamerInteractionEnergies, rotamers);

                // get actual energy
                Peptide thisPeptide = Rotamer.reconstitute(peptide, rotamers);
                new GaussianInputFile(thisPeptide).write(String.format("test_peptides/a_star_%04d.gjf", count));
                count++;
                List<Interaction> interactions = OPLScalculator.getInteractions(thisPeptide);
                Double[][] energyMatrix = AStarEnergyCalculator.getRotamerEnergyMatrix(thisPeptide, interactions);

                compareBackboneAtoms(peptide, thisPeptide);

                double thisBackboneEnergy = Interaction.getBackboneEnergy(energyMatrix);
                System.out.printf("     backbone   predicted %12.4f   actual %12.4f   diff %6.4f\n", backboneEnergy, thisBackboneEnergy, thisBackboneEnergy-backboneEnergy);

                double thisSingleEnergies = 0.0;
                for (Rotamer rotamer : rotamers)
                    {
                        // get the rotamer energy without the reference energy
                        Double singleEnergy = Interaction.getRotamerEnergy(rotamer, energyMatrix, false);
                        thisSingleEnergies += singleEnergy;
                    
                        double predictedSingleEnergy = rotamerSelfEnergies.get(rotamer);
                        System.out.printf("     %-25s   predicted %12.4f   actual %12.4f   diff %6.4f\n", rotamer.description, predictedSingleEnergy, singleEnergy, singleEnergy-predictedSingleEnergy);
                    }

                double thisInteractionEnergies = 0.0;
                for (int i=0; i < rotamers.size()-1; i++)
                    {
                        Rotamer rotamer1 = rotamers.get(i);
                        for (int j=i+1; j < rotamers.size(); j++)
                            {
                                Rotamer rotamer2 = rotamers.get(j);
                                RotamerPair pair = new RotamerPair(rotamer1, rotamer2);
                                double doubleEnergy = OPLScalculator.getInteractionEnergy(rotamer1, null, rotamer2, null);
                                thisInteractionEnergies += doubleEnergy;
                                double predictedDoubleEnergy = rotamerInteractionEnergies.get(pair);

                                System.out.printf("     %-50s   predicted %12.4f   actual %12.4f   diff %6.4f\n", rotamer1.description + " / " + rotamer2.description, predictedDoubleEnergy, doubleEnergy, doubleEnergy-predictedDoubleEnergy);
                            }
                    }
                
                double actualEnergy = thisBackboneEnergy + thisSingleEnergies + thisInteractionEnergies;

                // compare
                System.out.printf("predicted: %10.2f   actual: %10.2f   difference: %7.4f\n", predictedEnergy, actualEnergy, predictedEnergy-actualEnergy);
                double tempEnergy = 0.0;
                for (Interaction interaction : interactions)
                    tempEnergy += interaction.interactionEnergy;
                System.out.printf("total energy: %12.4f     diff %12.4f\n", tempEnergy, predictedEnergy-tempEnergy);
                System.out.print("\n\npress enter");
                Scanner scanner = new Scanner(System.in);
                scanner.nextLine();            
            }
    }
}
