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
 * Note that in the DEE classes, we have to consider the HN as part of the sidechain because proline doesn't have one.
 * In other words, we define the backbone as the set of atoms common to all conformations.  This is done by having the
 * Rotamer objects contain sidechain information only and the extra HN atom is added when considering rotamer pairs
 * in OPLScalculator.  Similarly, we include the extra HN atom in the sidechains in Interaction.java.<p>
 *
 * Energies are calculated over a truncated OPLS forcefield which include vdw interactions (radii scaled by 0.9) and
 * distance-dependent electrostatics (as defined in OPLScalculator and Interaction).
 * Reference energies are defined for each residue self-energy as the mean value of the residue self-energy (self-residue
 * + backbone-residue) over a set of unfolded and minimized peptide states.
 */
public class DEEenergyCalculator
{
    /** the self-energy of the sidechain and the interaction energy with the backbone, corrected by reference energy */
    public final Map<Rotamer,Double> rotamerSelfEnergies;

    /** the interaction energies between the two rotamers */
    public final Map<RotamerPair,Double> rotamerInteractionEnergies;

    /** the backbone self energy */
    public final double backboneEnergy;

    /** constructor */
    private DEEenergyCalculator(Map<Rotamer,Double> rotamerSelfEnergies,
                                Map<RotamerPair,Double> rotamerInteractionEnergies, double backboneEnergy)
    {
        this.rotamerSelfEnergies = rotamerSelfEnergies;
        this.rotamerInteractionEnergies = rotamerInteractionEnergies;
        this.backboneEnergy = backboneEnergy;
    }

    /**
     * Factory method to create a DEEenergyCalculator.  Calculates the single and pairwise energies of all the rotamers
     * in the RotamerSpace.  Energies on a truncated OPLS forcefield containing only charge and vdw interactions.
     * @param rotamerSpace to access the rotamer space and incompatible pairs
     * @return the energy calculations as contained in a DEEenergyCalculator object
     */
    public static DEEenergyCalculator analyze(RotamerSpace rotamerSpace)
    {
        // estimate the number of rotamers and rotamer pairs
        int totalRotamers = RotamerSpace.countRotamers(rotamerSpace.rotamerSpace);
        int estimated = ( ( totalRotamers * (totalRotamers - 1) ) / 2 ) - rotamerSpace.incompatiblePairs.size();

        // populate fields
        ConcurrentHashMap<Rotamer,Double> tempMap1 = new ConcurrentHashMap<Rotamer,Double>(totalRotamers, 0.75F, 16);
        ConcurrentHashMap<RotamerPair,Double> tempMap2 = new ConcurrentHashMap<RotamerPair,Double>(estimated, 0.75F, 16);
        AtomicDouble backboneEnergy = new AtomicDouble();
        
        // run energy jobs
        List<Future<Result>> futures = createEnergyJobs(rotamerSpace.rotamerSpace, rotamerSpace.peptide,
                                                        rotamerSpace.incompatiblePairs, tempMap1, tempMap2, backboneEnergy);
        GeneralThreadService.silentWaitForFutures(futures);

        // return the object
        return new DEEenergyCalculator(tempMap1, tempMap2, backboneEnergy.get());
    }

    /**
     * Debugging method that predicts the energy of a peptide constructed from a list of rotamers.
     * Will give erroneous results if there aren't enough rotamers to fill every position.  Hairpin positions don't count,
     * since they're defined as part of the backbone.
     * @param rotamerList the rotamers
     * @return the predicted energy in kcal
     */
    public double predictEnergy(RotamerSpace rotamerSpace, List<Rotamer> rotamerList)
    {
        // add the backbone energy
        double predictedEnergy = backboneEnergy;
        System.out.printf("backbone : %.4f\n", backboneEnergy);

        // add the energies of all the rotamers
        for (Rotamer r : rotamerList)
            {
                int index = rotamerSpace.rotamerSpace.get(r.sequenceIndex).indexOf(r);
                if ( ! rotamerSelfEnergies.containsKey(r) )
                    throw new IllegalArgumentException("missing energy for rotamer " + r.toString() + "  index = " + index);
                Double energy = rotamerSelfEnergies.get(r);
                predictedEnergy += energy;
                System.out.printf("%25s (%2d) : %12.4f\n", r.description, r.sequenceIndex, energy);
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
                        System.out.printf("%25s (%2d) / rotamer %25s (%2d) : %12.4f\n", rotamer1.description, rotamer1.sequenceIndex,
                                          rotamer2.description, rotamer2.sequenceIndex, energy);
                    }
            }

        return predictedEnergy;
    }

    /**
     * Creates all the RotamerEnergyJobs we need to get rotamer and rotamer pair energies
     * for an entire rotamer space.
     * @param rotamerSpace all the rotamers to analyze
     * @param startingPeptide contains the peptide backbone
     * @param incompatiblePairs all the incompatible pairs of rotamers we don't need interaction energies for
     * @param tempMap1 the map to concurrently update for rotamer single energies
     * @param tempMap2 the map to concurrently update for rotamer pair energies
     * @param backboneEnergy the backbone energy
     * @return the results of all the calculations (dummy objects in this case)
     */
    public static List<Future<Result>> createEnergyJobs(List<List<Rotamer>> rotamerSpace, Peptide startingPeptide, Set<RotamerPair> incompatiblePairs,
                                                        ConcurrentHashMap<Rotamer,Double> tempMap1, ConcurrentHashMap<RotamerPair,Double> tempMap2,
                                                        AtomicDouble backboneEnergy)
    {
        //System.out.println(new Date());
        //System.out.print("Preparing to create energy jobs...");

        // create all glycine peptide
        Peptide templatePeptide = startingPeptide;
        ProtoAminoAcid glycine = ProtoAminoAcidDatabase.getTemplate("standard_glycine");
        for (Residue residue : templatePeptide.sequence)
            {
                if ( residue.isHairpin )
                    continue;
                templatePeptide = SidechainMutator.mutateSidechain(templatePeptide, residue, glycine);
            }

        // create jobs
        //System.out.print("single rotamer jobs...");
        List<Future<Result>> futures = new ArrayList<>(); 
        
        // create rotamer jobs
        // we minimize the amount of work by going across each row of rotamers
        // that is, we take the first rotamer at each position and make a peptide
        // then we go down the list, skipping positions where there aren't enough rotamers
        // incompatibles should not be an issue because those interactions will be ignored
        int maxIndex = 0;
        for (List<Rotamer> list : rotamerSpace)
            {
                if ( list.size() - 1 > maxIndex )
                    maxIndex = list.size() -1;
            }

        for (int i=0; i <= maxIndex; i++)
            {
                List<Rotamer> theseRotamers = new LinkedList<>();
                for (List<Rotamer> list : rotamerSpace)    
                    {
                        if ( i < list.size() )
                            theseRotamers.add( list.get(i) );
                    }
                ImmutableList<Rotamer> theseRotamers2 = ImmutableList.copyOf(theseRotamers);
                RotamerEnergyJob job = new RotamerEnergyJob(templatePeptide, theseRotamers2, theseRotamers2, tempMap1, backboneEnergy);
                Future<Result> f = GeneralThreadService.submit(job);
                futures.add(f);
            }

        // create rotamer pair jobs
        // we save time by just computing the interactions between the rotamers without the backbone present
        //System.out.println(new Date());
        //System.out.println("rotamer pair jobs...");
        AStarEnergyCalculator.PairIterator iterator = new AStarEnergyCalculator.PairIterator(rotamerSpace, 1000);
        int count = 0;
        while (iterator.hasNext())
            {
                count++;
                LinkedListMultimap<Rotamer,Rotamer> thisBatch = iterator.next();
                RotamerPairEnergyJob job = new RotamerPairEnergyJob(templatePeptide, incompatiblePairs, thisBatch, tempMap2);
                Future<Result> f = GeneralThreadService.submit(job);
                futures.add(f);
                //System.out.printf("%d jobs created   \r", count);
            }
        //System.out.println("\nAll jobs created.");
        return futures;
    }

    /**
     * Calculates the self-energies of the specified peptides.
     */
    public static class RotamerEnergyJob implements WorkUnit
    {
        /** the starting peptide */
        public final Peptide startingPeptide;

        /** the rotamers to use to create the peptide we need to analyze */
        public final ImmutableList<Rotamer> rotamersToReconstitute;

        /** the rotamers to obtain the energies of */
        public final ImmutableList<Rotamer> rotamers;

        /** where to put the rotamer energies */
        public final ConcurrentHashMap<Rotamer,Double> map;

        /** where to put the backbone energy */
        public final AtomicDouble backboneEnergy;

        /**
         * Constructor.  We don't check that everything is consistent.  In particular, we don't check that rotamerToReconstitute
         * is compatible with rotamers and rotamerPairs.
         * @param startingPeptide the starting peptide, which might not have the correct rotamers
         * @param rotamersToReconstitute all the rotamers we need to use to create the peptide that will have the righht rotamers
         * @param rotamers the rotamers to compute the energies of
         * @param map the map to concurrently update with single rotamer energies
         * @param backboneEnergy the backbone energy
         */
        public RotamerEnergyJob(Peptide startingPeptide, ImmutableList<Rotamer> rotamersToReconstitute,
                                ImmutableList<Rotamer> rotamers, ConcurrentHashMap<Rotamer,Double> map, AtomicDouble backboneEnergy)
        {
            this.startingPeptide = startingPeptide;
            this.rotamersToReconstitute = rotamersToReconstitute;
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
            Peptide peptide = Rotamer.reconstitute(startingPeptide, rotamersToReconstitute);

            // analyze interactions
            List<Interaction> interactions = OPLScalculator.getInteractions(peptide);

            // initialize maps that will contain the answers
            Map<Rotamer,Double> rotamerMap = new LinkedHashMap<>();

            // analyze all interactions in this peptide
            Double[][] energyMatrix = Interaction.getRotamerEnergyMatrix(peptide, interactions);

            // get rotamer energies
            for (Rotamer rotamer : rotamers)
                {
                    Double singleEnergy = Interaction.getRotamerEnergy(rotamer, energyMatrix); 
                    rotamerMap.put(rotamer,singleEnergy);
                }

            // get backbone energy
            double thisBackboneEnergy = Interaction.getBackboneEnergy(energyMatrix); 

            // concurrently update results
            map.putAll(rotamerMap);
            backboneEnergy.compareAndSet(0.0, thisBackboneEnergy);

            // return result
            //return new RotamerEnergyResult(rotamerMap, rotamerPairMap, thisBackboneEnergy); 
            return null;
        }

        @Override
        public String toString()
        {
            return String.format("RotamerEnergyJob for %d rotamers", rotamers.size());
        }

        @Override
        public boolean equals(Object obj)
        {
            if ( obj == null )
                return false;
            if ( obj == this )
                return true;
            if ( !(obj instanceof RotamerEnergyJob) )
                return false;

            RotamerEnergyJob j = (RotamerEnergyJob)obj;
            if ( startingPeptide.equals(j.startingPeptide) &&
                 rotamersToReconstitute.equals(j.rotamersToReconstitute) &&
                 rotamers.equals(j.rotamers) )
                return true;
            return false;
        }

        @Override
        public int hashCode()
        {
            return Objects.hash(startingPeptide, rotamersToReconstitute, rotamers);
        }
    }

    /**
     * Calculates the interaction energies for the specified rotamers.
     */
    public static class RotamerPairEnergyJob implements WorkUnit
    {
        /** the peptide */
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
         * @return a dummy result; the energies will be updated concurrently through the concurrent map
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
                    
                    Atom extraAtom1 = null;
                    if ( rotamer1.description.indexOf("proline") == -1 )
                        extraAtom1 = peptide.sequence.get(rotamer1.sequenceIndex).HN;
                    
                    Atom extraAtom2 = null;
                    if ( rotamer2.description.indexOf("proline") == -1 )
                        extraAtom2 = peptide.sequence.get(rotamer2.sequenceIndex).HN;
                    
                    double energy = OPLScalculator.getInteractionEnergy(rotamer1, extraAtom1, rotamer2, extraAtom2);
                    tempMap.put(pair,energy);
                }
            map.putAll(tempMap);
            return null;
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
        int forbiddenIndex = (sequenceLength/2) - 1;
        List<String> stringSequence = ImmutableList.of("glycine", "standard_alanine", "standard_alanine", "glycine", "glycine", "glycine",
                                                       "standard_alanine", "glycine", "serine", "glycine", "arginine", "glycine");
        List<ProtoAminoAcid> protoAminoAcids = ProtoAminoAcidDatabase.getSpecificSequence(stringSequence);
        int tempJ = 0;
        for (int i=0; i < sequenceLength; i++)
            {
                if ( i == forbiddenIndex || i == forbiddenIndex+1 )
                    continue;
                Residue residue = peptide.sequence.get(i);
                ProtoAminoAcid protoAminoAcid = protoAminoAcids.get(tempJ);
                peptide = SidechainMutator.mutateSidechain(peptide, residue, protoAminoAcid);
                tempJ++; 
            }
        new GaussianInputFile(peptide).write("test_peptides/original.gjf");

        // create the rotamer space
        CatalystRotamerSpace catalystRotamerSpace = new CatalystRotamerSpace(peptide, true); 
        List<List<Rotamer>> rotamerSpace = catalystRotamerSpace.rotamerSpace;
        System.out.println("=== Rotamer Space: ===");
        for (int i=0; i < rotamerSpace.size(); i++)
            {
                System.out.println("Position " + i + ":");
                List<Rotamer> list = rotamerSpace.get(i);
                System.out.print("    ");
                for (Rotamer r : list)
                    System.out.print(r + "  ");
                System.out.println("\n");
            }

        // get energy of all rotamers and rotamer pairs
        DEEenergyCalculator calculator = DEEenergyCalculator.analyze(catalystRotamerSpace);

        // choose a random set of rotamers that don't include any incompatible pairs
        System.out.println("Drawing random rotamers...");
        
        Set<RotamerPair> incompatible = catalystRotamerSpace.incompatiblePairs;
        ThreadLocalRandom random = ThreadLocalRandom.current();
        for (int iteration=0; iteration < 150; iteration++)
            {
                // choose a random set of rotamers
                List<Rotamer> randomRotamers = null;
                System.out.println("\n\n\nRolling the dice... ");
                rolling:
                while (true)
                    {
                        randomRotamers = new LinkedList<>();
                        
                        // choose a random rotamer at each position
                        //hasProline = false;
                        for (List<Rotamer> list : rotamerSpace)
                            {
                                if ( list.size() == 0 )
                                    continue;
                                Rotamer randomRotamer = list.get( random.nextInt(list.size()) );
                                randomRotamers.add(randomRotamer);
                            }

                        // check for incompatible pairs
                        for (int i=0; i < randomRotamers.size(); i++)
                            {
                                Rotamer rotamer1 = randomRotamers.get(i);
                                for (int j=i+1; j < randomRotamers.size(); j++)
                                    {
                                        Rotamer rotamer2 = randomRotamers.get(j);
                                        RotamerPair pair = new RotamerPair(rotamer1, rotamer2);
                                        if ( incompatible.contains(pair) )
                                            {
                                                System.out.print("rejected ");
                                                continue rolling;
                                            }
                                    }
                            }

                        System.out.println("=== rotamer list ===");
                        for (Rotamer r : randomRotamers)
                            System.out.println(r);
                        System.out.println(">>> accepted");
                        break;
                    }

                // reconstitute the peptide
                Peptide reconstitutedPeptide = Rotamer.reconstitute(peptide, randomRotamers);

                GaussianInputFile f = new GaussianInputFile(reconstitutedPeptide);
                f.write("test_peptides/reconstituted.gjf");

                // check the actual energy
                List<Interaction> interactions = new ArrayList<>(OPLScalculator.getInteractions(reconstitutedPeptide));
                
                System.out.println("\ntop interactions:\n");
                Collections.sort(interactions);
                for (int i=0; i < 20; i++)
                    System.out.println(interactions.get(i).toString(reconstitutedPeptide));
                System.out.println();

                Double[][] energyMatrix = Interaction.getRotamerEnergyMatrix(reconstitutedPeptide, interactions);
                System.out.println("\nbreakdown of actual energy:");
                double actualBackboneEnergy = energyMatrix[energyMatrix.length-1][energyMatrix.length-1];
                for (int i=0; i <= reconstitutedPeptide.sequence.size(); i++)
                    {
                        String string1 = "";
                        if ( i < reconstitutedPeptide.sequence.size() )
                            string1 = reconstitutedPeptide.sequence.get(i).description;
                        else
                            string1 = "backbone";
                        for (int j=0; j <= reconstitutedPeptide.sequence.size(); j++)
                            {
                                String string2 = "";
                                if ( j < reconstitutedPeptide.sequence.size() )
                                    string2 = reconstitutedPeptide.sequence.get(j).description;
                                else
                                    string2 = "backbone";
                                System.out.printf("%s(%d) - %s(%d) : %.4f\n", string1, i, string2, j, energyMatrix[i][j]);
                            }
                    }
                double totalEnergy = 0.0;
                for (Interaction i : interactions)
                    totalEnergy += i.interactionEnergy;
                
                // add reference energy correction
                for (Rotamer r : randomRotamers)
                    totalEnergy -= Interaction.getReferenceEnergy(r);

                double predictedEnergy = calculator.predictEnergy(catalystRotamerSpace, randomRotamers);
                double difference = (predictedEnergy - totalEnergy);
                System.out.printf("predicted total energy: %.4f\n", predictedEnergy);
                System.out.printf("actual total energy: %.4f\n\n", totalEnergy);
                if ( Math.abs(difference) > 0.1 )
                    System.out.printf("*** ");
                System.out.printf("difference: %.4f\n", difference);

                //if ( Math.abs(difference) > 0.1 )
                if ( true )
                    {
                        System.out.println("================");
                        for (Rotamer r : randomRotamers)
                            {
                                double thisPredictedEnergy = calculator.rotamerSelfEnergies.get(r);
                                double thisActualEnergy = energyMatrix[r.sequenceIndex][r.sequenceIndex] + energyMatrix[r.sequenceIndex][energyMatrix.length-1] -
                                                          Interaction.getReferenceEnergy(r);
                                double thisDifference = thisActualEnergy - thisPredictedEnergy;
                                System.out.printf("%-3d %-40s  pred: %12.4f   actual: %12.4f   diff: %12.4f\n",
                                                 r.sequenceIndex, r.description, thisPredictedEnergy, thisActualEnergy, thisDifference);
                            
                                if ( Math.abs(difference) > 0.1 )
                                   {
                                        reconstitutedPeptide = Rotamer.reconstitute(peptide, ImmutableList.of(r));
                                        List<Interaction> extraInteractions = new ArrayList<>(OPLScalculator.getInteractions(reconstitutedPeptide));
                                        Double[][] extraMatrix = Interaction.getRotamerEnergyMatrix(reconstitutedPeptide, extraInteractions);
                                        int thisIndex = r.sequenceIndex;
                                        System.out.println(extraMatrix[thisIndex][thisIndex]+extraMatrix[thisIndex][extraMatrix.length-1]);
                                    }
                            }
                        for (int i=0; i < randomRotamers.size(); i++)
                            {
                                Rotamer rotamer1 = randomRotamers.get(i);
                                for (int j=i+1; j < randomRotamers.size(); j++)
                                    {
                                        Rotamer rotamer2 = randomRotamers.get(j);
                                        RotamerPair pair = new RotamerPair(rotamer1, rotamer2);
                                        double thisPredictedEnergy = calculator.rotamerInteractionEnergies.get(pair);
                                        double thisActualEnergy = energyMatrix[rotamer1.sequenceIndex][rotamer2.sequenceIndex];
                                        double thisDifference = thisActualEnergy - thisPredictedEnergy;
                                        System.out.printf("%-2d / %-2d %-40s / %-40s  pred: %12.4f   actual: %12.4f   diff: %12.4f\n",
                                                 rotamer1.sequenceIndex, rotamer2.sequenceIndex,
                                                 rotamer1.description, rotamer2.description,
                                                 thisPredictedEnergy, thisActualEnergy, thisDifference);
                                    }
                            }
                        System.out.println("================");
                    }

                // wait for user input
                //if ( Math.abs(difference) > 0.1 )
                if ( true )
                    {
                        System.out.print("\n\npress enter");
                        Scanner scanner = new Scanner(System.in);
                        scanner.nextLine();
                    }
            }
    }
}
