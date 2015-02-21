import java.util.*;
import java.util.concurrent.*;
import java.util.concurrent.atomic.*;
import java.io.*;
import com.google.common.collect.*;

/**
 * This class calculates reference energies for amino acids in a beta sheet confromation. 
 * Reference energies are collected for both OPLS (calculated internally) and AMOEBA (calculated
 * externally by TINKER).  These energies do not include solvation energies.
 * 
 * The procedure is:
 * 
 * 1. Draw random poly-gly peptides in the beta sheet conformation using BetaSheetGenerator
 *
 * 2. Perform mutations to generate random beta sheet peptides using SidechainMutator
 *
 * 3. Perform a Monte Carlo minimization on each randomly generated peptide.  These minimizations
 *    are performed with OPLS in the gas phase.
 *
 * 4. Perform approximate solvation on the structures from step (3) to select the best structures.
 *
 * 5. Minimize the best structures with AMOEBA in the gas phase.  Add a single point solvation correction
 *    with GK at the end.  Take the lowest energy structure for each peptide and put it in a pile.
 *
 * 6. For all the structures in the pile, compute the gas phase AMOEBA energy breakdown and use that for the
 *    AMOEBA reference energies.  Also compute the OPLS interactions and use those for the OPLS energies.
 */
public class ReferenceEnergyCalculator
{
    /** The number of reference peptides to minimize. */
    public static final int NUMBER_OF_REFERENCE_PEPTIDES = 1;

    /** The number of poses per reference peptide to minimize. */
    public static final int NUMBER_OF_STRUCTURES_TO_MINIMIZE = 500;

    /** not instantiable */
    private ReferenceEnergyCalculator()
    {
        throw new IllegalArgumentException("cannot be instantiated!");
    }

    /**
     * This method will perform mutations to a list of beta sheet peptides to randomly generate a list of peptides. 
     * We mutate to hairpins containing exactly one his, arg, and ser on the same face.  The resulting peptides may
     * clash.  Iterates through the input list as many times as needed to make the target number of results.
     * It will also save these peptides to disk.
     * @param betaSheets the initial beta sheets to be mutated
     * @param targetNumberOfResults the number of peptides to generate
     * @return a list of random sequence peptides still in the beta sheet secondary structure
     */
    public static List<Peptide> generateRandomPeptides(List<Peptide> betaSheets, int targetNumberOfResults)
    {
        // verify invariatns
        if ( betaSheets == null || betaSheets.size() < 1 )
            throw new IllegalArgumentException("check beta sheets");
        if ( targetNumberOfResults < 1 )
            throw new IllegalArgumentException("check target number of results");

        // get templates
        List<ProtoAminoAcid> mutationOutcomes = CatalystRotamerSpace.MUTATION_OUTCOMES;
        ProtoAminoAcid arginine = ProtoAminoAcidDatabase.getTemplate("arginine");
        ProtoAminoAcid histidineHD = ProtoAminoAcidDatabase.getTemplate("histidine_hd");
        ProtoAminoAcid histidineHE = ProtoAminoAcidDatabase.getTemplate("histidine_he");
        ProtoAminoAcid serine = ProtoAminoAcidDatabase.getTemplate("serine");

        // perform mutations
        List<Peptide> randomPeptides = new ArrayList<>(targetNumberOfResults);
        int peptideCount = 0;
        while (true)
            {
                // get the current template
                if ( randomPeptides.size() >= targetNumberOfResults )
                    break;
                Peptide p = betaSheets.get(peptideCount);
                peptideCount++;
                if ( peptideCount > betaSheets.size() - 1 )
                    peptideCount = 0;

                // determine where mutations can be made
                LinkedList<Integer> allowedPositions = new LinkedList<>();
                for (int i=0; i < p.sequence.size(); i++)
                    {
                        if ( p.sequence.get(i).isHairpin )
                            continue;
                        allowedPositions.add(i);
                    }

                // place arginine
                Peptide peptide = p;
                Collections.shuffle(allowedPositions);
                int randomIndex = allowedPositions.remove();
                Residue residue = peptide.sequence.get(randomIndex);
                peptide = SidechainMutator.mutateSidechain(peptide, residue, arginine);
                
                // place histidine and serine on the same side
                boolean isUp = RotamerFactory.isUp(p.sequence.size(), randomIndex);
                List<Integer> toBeRemoved = new ArrayList<>(2);
                for (Integer i : allowedPositions)
                    {
                        boolean thisIsUp = RotamerFactory.isUp(p.sequence.size(), i);
                        if ( isUp == thisIsUp )
                            {
                                toBeRemoved.add(i);
                                residue = peptide.sequence.get(i);
                                
                                // pick a random histidine
                                ProtoAminoAcid randomHistidine = histidineHD;
                                if ( ThreadLocalRandom.current().nextDouble() < 0.50 )
                                    randomHistidine = histidineHE;
                                peptide = SidechainMutator.mutateSidechain(peptide, residue, randomHistidine);
                                break;
                            }
                    }
                 for (Integer i : allowedPositions)
                    {
                        boolean thisIsUp = RotamerFactory.isUp(p.sequence.size(), i);
                        if ( isUp == thisIsUp )
                            {
                                toBeRemoved.add(i);
                                residue = peptide.sequence.get(i);
                                peptide = SidechainMutator.mutateSidechain(peptide, residue, serine);
                                break;
                            }
                    }
                if ( toBeRemoved.size() != 2 )
                    throw new IllegalArgumentException("unreachable");
                allowedPositions.removeAll(toBeRemoved);

                // mutate all the other positions at random
                for (Integer i : allowedPositions)
                    {
                        randomIndex = ThreadLocalRandom.current().nextInt(mutationOutcomes.size());
                        ProtoAminoAcid randomTemplate = mutationOutcomes.get(randomIndex);
                        residue = peptide.sequence.get(i);
                        peptide = SidechainMutator.mutateSidechain(peptide, residue, randomTemplate);
                    }

                // add to result
                randomPeptides.add(peptide);
            }
        return randomPeptides;
    }

    /** 
     * Sorts the energies by residue from a bunch of peptides by amino acid description.
     */
    public static Map<String,List<Double>> getEnergies(List<Peptide> peptides)
    {
        Map<String,List<Double>> returnMap = new HashMap<>();
        for (Peptide p : peptides)
            {
                int sequenceLength = p.sequence.size();
                List<Double> energiesByResidue = p.energyBreakdown.energyByResidue;
                for (int i=0; i < sequenceLength; i++)
                    {
                        Residue residue = p.sequence.get(i);
                        if ( residue.isHairpin )
                            continue;
                        String description = residue.description;
                        double energy = energiesByResidue.get(i);
                        List<Double> list = returnMap.get(description);
                        if ( list == null )
                            {
                                list = new ArrayList<Double>();
                                returnMap.put(description, list);
                            }
                        list.add(energy);
                    }
            }
        return returnMap;
    }

//    /** 
//    * A method to make the individual calls needed to run the process of finding reference energies
//    * @return a map from amino acids to reference energies
//    */
//    public static Map<AminoAcid, Double> calculateReferenceEnergies()
//    {
//        // Draw random beta sheet backbones
//        // Generate 5000 backbones and pick 500 lowest energy ones
//        
//        List<Peptide> sheets = BetaSheetGenerator.generateSheets(5, 1000, 5000, .01);
//        Collections.sort(sheets);
//        List<Peptide> lowEnergyBackbones = sheets;
//        
//        // Pick lowest energy structures
//        for (int i = 0; i < 500; i++)
//            lowEnergyBackbones.add(sheets.get(i));
//        
//        // Perform mutations
//        List<Peptide> randomPeptides = generateRandomPeptides(lowEnergyBackbones);
//
//        // Launch Monte Carlo jobs
//        List<Peptide> minimizedRandomPeptides = new LinkedList<>();
//        for (Peptide p : randomPeptides)
//        {
//            List<Peptide> monteCarloOutput = runMonteCarlo(p);
//            // Minimize with AMOEBA
//            List<Peptide> lowestEnergyMinimizedPeptides = minimize(monteCarloOutput);
//            for (Peptide p2 : lowestEnergyMinimizedPeptides)
//                minimizedRandomPeptides.add(p2);
//        }
//
//        // Breakdown energy for each peptide and average energies
//        Map<AminoAcid, Double> referenceEnergies = averageEnergies(minimizedRandomPeptides);
//        return referenceEnergies;
//    }
//
    /**
     * Performs the reference energy generation procedure in parallel on one node.
     */
    public static void main(String[] args)
    {
        // create some beta sheets
        DatabaseLoader.go();
        List<Peptide> sheets = BetaSheetGenerator.generateSheets(5,     // arm length 
                                                                 10,    // max iterations
                                                                 100,   // max results
                                                                 0.01); // cooling rate
        System.out.printf("%d beta sheets have been generated.\n", sheets.size());

        // randomly mutate
        List<Peptide> initialRandomPeptides = generateRandomPeptides(sheets, NUMBER_OF_REFERENCE_PEPTIDES);
        //for (Peptide p : initialRandomPeptides)
        //    System.out.println(p.name);

        // minimize with OPLS
        List<Peptide> minimizedRandomPeptides = TinkerJob.minimize(initialRandomPeptides, Forcefield.OPLS, 3000, false, false, false, false, false);
        System.out.printf("\n%d initial structues generated.\n", minimizedRandomPeptides.size());
        //Peptide.writeGJFs(minimizedRandomPeptides, "test_peptides/random_", 3, 10);

        // check for backbone clashes
        List<Peptide> startingPeptides = new ArrayList<>();
        for (int i=0; i < minimizedRandomPeptides.size(); i++)
            {
                System.out.printf("Checking clashes %d of %d...\r", i+1, startingPeptides.size());
                Peptide p = minimizedRandomPeptides.get(i);
                List<Pair<Integer,Integer>> backbonePairs = p.getBackbonePairs();
                if (!p.hasBackboneClash(backbonePairs))
                    startingPeptides.add(p);
           }
        System.out.printf("\nDone.  %d peptides passed the clash check.\n", startingPeptides.size());

        // run the monte carlo jobs
        // --->>>>>>>> add checkpoint filenames
        List<Future<Result>> futures = new ArrayList<>(startingPeptides.size());
        for (Peptide p : startingPeptides)
            {
                FixedSequenceMonteCarloJob job = new FixedSequenceMonteCarloJob(p, 0.001, 1000, NUMBER_OF_STRUCTURES_TO_MINIMIZE, 4, null);
                Future<Result> f = GeneralThreadService.submit(job);
            }
        GeneralThreadService.silentWaitForFutures(futures);

        // get the results
        List<List<Peptide>> bestPoses = new ArrayList<>(futures.size()); // outer list is indexed by starting peptide, inner list by pose
        for (Future<Result> f : futures)
            {
                MonteCarloResult result = null;
                try { result = (MonteCarloResult)f.get(); }
                catch (Exception e) { e.printStackTrace(); continue; }
                bestPoses.add(result.bestPeptides);
            }

        // minimize all poses with AMOEBA with single point GK solvation
        List<List<Peptide>> bestPoses2 = new ArrayList<>(futures.size());
        for (List<Peptide> list : bestPoses)
            {
                List<Peptide> list2 = TinkerJob.minimize(list, Forcefield.AMOEBA, 2000, false, false, true, true, false); 
                bestPoses2.add(list2);
            }

        // take the best pose for each peptide and perform a gas phase AMOEBA energy breakdown
        List<Peptide> bestPoses3 = new ArrayList<>(bestPoses2.size());
        for (List<Peptide> list : bestPoses2)
            {
                if ( list.size() > 0 )
                    bestPoses3.add(list.get(0));
            }
        List<Peptide> bestPoses4 = TinkerJob.analyze(bestPoses3, Forcefield.AMOEBA);
        
        // get AMOEBA reference energies
        Map<String,List<Double>> AMOEBAreferenceEnergies = getEnergies(bestPoses4);
        for (String description : AMOEBAreferenceEnergies.keySet())
            {
                List<Double> energies = AMOEBAreferenceEnergies.get(description);
                System.out.printf("%s (%d energies)\n", description, energies.size());
                System.out.println(energies);
            }

        // take the best poses for each peptide and perform a gas phase OPLS interaction calculation

        // get OPLS reference energies
    }
}
