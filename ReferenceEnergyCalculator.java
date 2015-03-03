import java.util.*;
import java.util.concurrent.*;
import java.util.concurrent.atomic.*;
import java.io.*;
import com.google.common.collect.*;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;

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
    public static final int NUMBER_OF_REFERENCE_PEPTIDES = 64;

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
        // verify invariants
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
                int charge = PeptideChargeCalculator.getCharge(peptide);
                if ( charge == 0 || charge == 1 )
                    randomPeptides.add(peptide);
            }
        return randomPeptides;
    }

    /** 
     * Sorts the energies by residue from a bunch of peptides by amino acid description.
     */
    public static Map<String,List<Double>> getAMOEBAEnergies(List<Peptide> peptides)
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

    /**
     * Computes one-center OPLS energies for DEE.
     */
    public static Map<String,List<Double>> getOPLSEnergies(List<Peptide> peptides)
    {
        Map<String,List<Double>> returnMap = new HashMap<>();
        for (Peptide p : peptides)
            {
                int sequenceLength = p.sequence.size();

                // compute OPLS energies
                List<Interaction> interactions = new ArrayList<>(OPLScalculator.getInteractions(p));
                Double[][] energyMatrix = Interaction.getRotamerEnergyMatrix(p, interactions, true);
                List<Double> energiesByResidue = new ArrayList<>(sequenceLength);
                for (int i=0; i < sequenceLength; i++)
                    {
                        double thisEnergy = energyMatrix[i][i] + energyMatrix[i][energyMatrix.length-1]; // rotamer self energy + rotamer/backbone energy
                        energiesByResidue.add(thisEnergy);
                    }

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

    /**
     * Performs the reference energy generation procedure in parallel on one node.
     */
    public static void main(String[] args)
    {
        // create some beta sheets
        DatabaseLoader.go();
        List<Peptide> sheets = BetaSheetGenerator.generateSheets(5,     // arm length 
                                                                 10,   // max iterations
                                                                 1000,  // max results
                                                                 0.05); // cooling rate
        System.out.printf("%d beta sheets have been generated.\n", sheets.size());

        // randomly mutate
        List<Peptide> initialRandomPeptides = generateRandomPeptides(sheets, NUMBER_OF_REFERENCE_PEPTIDES);
        //for (Peptide p : initialRandomPeptides)
        //    System.out.println(p.name);

        // minimize with OPLS
        System.out.println("Minimizing intial peptides...");
        List<Peptide> minimizedRandomPeptides = TinkerJob.minimize(initialRandomPeptides, Forcefield.OPLS, 3000, false, false, false, false, false);
        System.out.printf("\n%d initial structues generated.\n", minimizedRandomPeptides.size());
        //Peptide.writeGJFs(minimizedRandomPeptides, "test_peptides/random_", 3, 10);

        // check for backbone clashes
        List<Peptide> startingPeptides = new ArrayList<>();
        for (int i=0; i < minimizedRandomPeptides.size(); i++)
            {
                System.out.printf("Checking clashes %d of %d...\r", i+1, minimizedRandomPeptides.size());
                Peptide p = minimizedRandomPeptides.get(i);
                List<Pair<Integer,Integer>> backbonePairs = p.getBackbonePairs();
                if (!p.hasBackboneClash(backbonePairs))
                    startingPeptides.add(p);
           }
        if ( startingPeptides.size() == 0 )
            throw new IllegalArgumentException("no peptides to start with");
        System.out.printf("\nDone.  %d peptides passed the clash check.\n", startingPeptides.size());

        // rotamer pack all of the starting structures
        // make a list of initial rotamers of the starting structure and use the best structure as the Monte Carlo seed
        List<List<Peptide>> startingPeptides2 = new ArrayList<>(startingPeptides.size());
        for (Peptide p : startingPeptides)
            {
                // perform A* iteration
                try
                    {
                        FixedSequenceRotamerSpace fixedSequenceRotamerSpace = new FixedSequenceRotamerSpace(p, false);
                        AStarEnergyCalculator calculator = AStarEnergyCalculator.analyze(fixedSequenceRotamerSpace);
                        RotamerIterator iterator = new RotamerIterator(fixedSequenceRotamerSpace.rotamerSpace, calculator.rotamerSelfEnergies, calculator.rotamerInteractionEnergies, 100);
                        List<RotamerIterator.Node> solutions = iterator.iterate();
                        solutions = iterator.iterate();
                        if ( solutions.size() == 0 )
                                continue;

                        // reconstitute
                        List<Peptide> results = new ArrayList<>(100);
                        for (RotamerIterator.Node node : solutions)
                            {
                                if (results.size() >= 100)
                                    break;
                                results.add(Rotamer.reconstitute(p, node.rotamers));
                           }
                        
                        // minimize peptides in parallel
                        List<Peptide> results2 = TinkerJob.minimize(results, Forcefield.OPLS, 2000, false, false, true, true, false);
                        startingPeptides2.add(results2);
                        double bestEnergy = results2.get(0).energyBreakdown.totalEnergy;
                        System.out.printf("Initial rotamer packing done for peptide %d of %d.  %d structures added (bestE = %.2f).\n", startingPeptides2.size(), startingPeptides.size(), results2.size(), bestEnergy);
                    }
                catch (Exception e)
                    {
                        e.printStackTrace();
                        continue;
                    }
                /*// use this to skip the initial rotamer packing
                List<Peptide> results2 = new ArrayList<>();
                results2.add(p);
                startingPeptides2.add(results2);
                */
            }

        // run the monte carlo jobs
        List<Future<Result>> futures = new ArrayList<>(startingPeptides.size());
        int peptideCount = 0;
        for (List<Peptide> list : startingPeptides2)
            {
                Peptide p = list.get(0);
                String filename = String.format("checkpoints/fsmcjob_%05d.chk", peptideCount);
                peptideCount++;
                FixedSequenceMonteCarloJob job = new FixedSequenceMonteCarloJob(p, 0.001, 1000, NUMBER_OF_STRUCTURES_TO_MINIMIZE, 4, filename, list);
                Future<Result> f = GeneralThreadService.submit(job);
                futures.add(f);
            }
        GeneralThreadService.silentWaitForFutures(futures);

        // get the results
        System.out.println("All MC jobs are complete.");
        List<List<Peptide>> bestPoses = new ArrayList<>(futures.size()); // outer list is indexed by starting peptide, inner list by pose
        for (Future<Result> f : futures)
            {
                MonteCarloResult result = null;
                try { result = (MonteCarloResult)f.get(); }
                catch (Exception e) { e.printStackTrace(); continue; }
                bestPoses.add(result.bestPeptides);
            }

        // minimize all poses with AMOEBA with approximate single point solvation
        System.out.println("Minimizing all poses with AMOEBA...");
        List<List<Peptide>> bestPoses2 = new ArrayList<>(futures.size());
        Map<Future<Result>,Integer> futureMap = new HashMap<>();
        for (int i=0; i < bestPoses.size(); i++)
            {
                List<Peptide> list = bestPoses.get(i);
                for (Peptide p : list)
                    {
                        TinkerJob job = new TinkerJob(p, Forcefield.AMOEBA, 2000, false, false, true, false, true);
                        Future<Result> f = GeneralThreadService.submit(job);
                        futureMap.put(f, i);
                    }
                bestPoses2.add(new ArrayList<Peptide>());
            }
        futures = new ArrayList<>(futureMap.keySet());
        GeneralThreadService.waitForFutures(futures);
        for (Future<Result> f : futureMap.keySet())
            {
                Integer index = futureMap.get(f);
                TinkerJob.TinkerResult result = null;
                try
                    {
                        result = (TinkerJob.TinkerResult)f.get();
                    }
                catch (Exception e)
                    {
                        e.printStackTrace();
                        continue;
                    }
                List<Peptide> list = bestPoses2.get(index);
                list.add(result.minimizedPeptide);
            }

        for (List<Peptide> list : bestPoses2)
            Collections.sort(list);

        // take the best pose for each peptide
        List<Peptide> bestPoses3 = new ArrayList<>(bestPoses2.size());
        for (List<Peptide> list : bestPoses2)
            {
                if ( list.size() > 0 )
                    bestPoses3.add(list.get(0));
                if ( list.size() > 1 )
                    bestPoses3.add(list.get(1));
            }
        Peptide.writeGJFs(bestPoses3, "test_peptides/best_", 3, 1000);
        Peptide.writeCHKs(bestPoses3, "test_peptides/best_", 3, 1000);

        // get AMOEBA reference energies
        System.out.println("Getting AMOEBA reference energies...");
        Map<String,List<Double>> AMOEBAreferenceEnergies = getAMOEBAEnergies(bestPoses3);
        String AMOEBAstring = "";
        for (String description : AMOEBAreferenceEnergies.keySet())
            {
                List<Double> energies = AMOEBAreferenceEnergies.get(description);
                DescriptiveStatistics stats = new DescriptiveStatistics();
                for (Double d : energies)
                    stats.addValue(d);
                double stderr = stats.getStandardDeviation() / Math.sqrt(stats.getN());
                String thisSummary = String.format("%30s : mean = %8.2f   stdev = %8.2f   err = %8.2f   n = %4d\n", description,
                                                   stats.getMean(), stats.getStandardDeviation(), stderr, stats.getN());
                System.out.print(thisSummary);
                AMOEBAstring += String.format("%30s %.2f\n", description, stats.getMean());
            }
        InputFileFormat.writeStringToDisk(AMOEBAstring, "AMOEBA.txt");

        // take the best poses for each peptide and perform a gas phase OPLS interaction calculation
        System.out.println("Getting OPLS reference energies...");
        Map<String,List<Double>> OPLSreferenceEnergies = getOPLSEnergies(bestPoses3);
        String OPLSstring = "";
        for (String description : OPLSreferenceEnergies.keySet())
            {
                List<Double> energies = OPLSreferenceEnergies.get(description);
                DescriptiveStatistics stats = new DescriptiveStatistics();
                for (Double d : energies)
                    stats.addValue(d);
                double stderr = stats.getStandardDeviation() / Math.sqrt(stats.getN());
                String thisSummary = String.format("%30s : mean = %8.2f   stdev = %8.2f   err = %8.2f   n = %4d\n", description,
                                                   stats.getMean(), stats.getStandardDeviation(), stderr, stats.getN());
                System.out.print(thisSummary);
                OPLSstring += String.format("%30s %.2f\n", description, stats.getMean());
            }
        InputFileFormat.writeStringToDisk(OPLSstring, "OPLS.txt");
    }
}
