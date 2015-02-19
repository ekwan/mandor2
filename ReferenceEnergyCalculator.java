import java.util.*;
import java.util.concurrent.*;
import java.util.concurrent.atomic.*;
import java.io.*;
import com.google.common.collect.*;

/**
 * This class calculates reference energies for amino acids in a beta sheet confromation. 
 * 
 * The procedure is:
 * 
 * 1. Draw random poly-gly peptides in the beta sheet conformation using BetaSheetGenerator
 *
 * 2. Perform mutations to generate random beta sheet peptides using SidechainMutator
 *
 * 3. Perform a Monte Carlo minimization on each randomly generated peptide
 *
 * 4. For lowest energy structures, minimize using AMOEBA and then perform an energy breakdown
 *
 * 5. Average across all amino acid occurences, for reference energy for each amino acid
 */

public class ReferenceEnergyCalculator
{
    /** Define constants */

    /** The number of poly-gly beta sheets to generate as a starting point */
    public static final int NUMBER_BETA_SHEETS = 500;

    /** The number of lowest energy structures to keep for each starting random peptide following fixed sequence Monte Carlo */
    public static final int NUMBER_LOWEST_ENERGY_STRUCTURES_TO_KEEP = 10;

    /** The number of structures to minimize using AMOEBA following the fixed sequence Monte Carlo */
    public static final int NUMBER_STRUCTURES_TO_MINIMIZE = 1000;

    /** not instantiable */
    private ReferenceEnergyCalculator()
    {
        throw new IllegalArgumentException("cannot be instantiated!");
    }

    /**
    * This method will perform mutations to a list of beta sheet peptides to randomly generate a list of peptides. 
    * This method does not permit mutations to cysteine, proline, and methionine. 
    * It will also save these peptides to disk.
    * @param betaSheets the initial beta sheets to be mutated
    * @return a list of random sequence peptides still in the beta sheet secondary structure
    */
    public static List<Peptide> generateRandomPeptides(List<Peptide> betaSheets)
    {
        // Create list of possible mutations
        List<ProtoAminoAcid> mutationOutcomes = new ArrayList<>();
        for (AminoAcid a : ProtoAminoAcidDatabase.KEYS)
            {
                // reject unusual amino acids
                if ( a == AminoAcid.DPRO || a == AminoAcid.LPRO || a == AminoAcid.TS ||
                     a == AminoAcid.LYS || a == AminoAcid.CYS  || a == AminoAcid.MET     )
                    continue;
                int index = ProtoAminoAcidDatabase.KEYS.indexOf(a);
                List<ProtoAminoAcid> VALUES = ProtoAminoAcidDatabase.VALUES.get(index);
                for (ProtoAminoAcid paa : VALUES)
                    {
                        // reject transition states
                        // decide upon this
                        if (paa.residue.description.toLowerCase().indexOf("hairpin") > -1 )
                            continue;
                        mutationOutcomes.add(paa);
                    }
            }

        // Generate random peptides
        List<Peptide> randomPeptides = new LinkedList<>();
        for (Peptide p : betaSheets)
        {
            Peptide randomPeptide = p;
            for (Residue r : p.sequence)
            {
                if (r.description.toLowerCase().indexOf("hairpin") > -1)
                    continue;

                // Pick random mutation target
                Collections.shuffle(mutationOutcomes);
                ProtoAminoAcid targetPAA = mutationOutcomes.get(0);
                randomPeptide = SidechainMutator.mutateSidechain(randomPeptide, r, targetPAA);
            }

            randomPeptides.add(randomPeptide);
        }
        
        return randomPeptides;
    }

    
    /** 
    * This method will run the fixed sequence Monte Carlo process.
    * @param randomPeptide a randomly generated peptide from performMutations
    * @return a list of size NUMBER_STRUCTURES_TO_MINIMIZE that will contain the lowest energy peptides from a fixed sequence simulation
    */
    public static List<Peptide> runMonteCarlo(Peptide randomPeptide)
    {
        FixedSequenceMonteCarloJob fixedSequenceMonteCarloJob = new FixedSequenceMonteCarloJob(randomPeptide, .01, 1000, NUMBER_STRUCTURES_TO_MINIMIZE); 
        MonteCarloResult monteCarloResult = fixedSequenceMonteCarloJob.call();
        return monteCarloResult.bestPeptides;
    }

    /**
    * A method for performing AMOEBA minimizations that also removes duplicates and returns the lowest energy structures only
    * @param monteCarloPeptides the peptides generated from the Monte Carlo process which includes NUMBER_STRUCTURES_TO_MINIMIZE peptides
    * @return the NUMBER_LOWEST_ENERGY_STRUCTURES_TO_KEEP peptides that will be used to find the average amino acid reference energy 
    */
    public static List<Peptide> minimize(List<Peptide> monteCarloPeptides)
    {
        // Create hashamp to avoid duplicates
        Map<PeptideFingerprint, Peptide> peptidesWithoutDuplicates = new HashMap<>();   

        // Only add low energy peptides
        List<Peptide> lowEnergyPeptides = new ArrayList<>();
        for (Peptide p : monteCarloPeptides)
        {
            // Minimize each peptide on AMOEBA forcefield
            // Use Tinker approximate solvation for analyze 
            TinkerJob job = new TinkerJob(p, Forcefield.AMOEBA, 2000, false, true, true, false, true);
            TinkerJob.TinkerResult result = job.call();
            Peptide newPeptide = result.minimizedPeptide; 
            
            PeptideFingerprint fingerprint = new PeptideFingerprint(newPeptide, newPeptide.energyBreakdown.totalEnergy);
            
            // Verify that we are adding new peptide and 
            if (!peptidesWithoutDuplicates.containsKey(fingerprint))
            {
                if (lowEnergyPeptides.size() < NUMBER_LOWEST_ENERGY_STRUCTURES_TO_KEEP)
                {
                    lowEnergyPeptides.add(newPeptide);
                    Collections.sort(lowEnergyPeptides);
                    peptidesWithoutDuplicates.put(fingerprint, newPeptide);
                }
                else
                {
                    // if highest energy member in list is lower in energy than current peptide then do not add it to the list
                    if (lowEnergyPeptides.get(lowEnergyPeptides.size()-1).energyBreakdown.totalEnergy > newPeptide.energyBreakdown.totalEnergy)
                    {
                        peptidesWithoutDuplicates.put(fingerprint, newPeptide);
                    }
                    else
                    {
                        lowEnergyPeptides.add(newPeptide);
                        Collections.sort(lowEnergyPeptides);
                        peptidesWithoutDuplicates.put(fingerprint, newPeptide);
                    }
                }
            }
        }

        return lowEnergyPeptides;
    }


    /** 
    * This method will calculate the average energy of each amino acid by taking an energy breakdown of each of the provided structures
    * @param structures a random set of low-energy structures 
    * @return a map from amino acids to average energies
    */
    public static Map<AminoAcid, Double> averageEnergies(List<Peptide> structures)
    {
        Map<AminoAcid, List<Double>> allEnergiesByAminoAcid = new HashMap<>();
        for (Peptide p : structures)
        {
            for (Residue r : p.sequence)
            {
                if (allEnergiesByAminoAcid.containsKey(r.aminoAcid))
                {
                    List<Double> energies = allEnergiesByAminoAcid.get(r.aminoAcid);
                    energies.add(p.energyBreakdown.energyByResidue.get(p.sequence.indexOf(r)));
                    allEnergiesByAminoAcid.put(r.aminoAcid, energies);
                }
                else
                {
                    List<Double> energies = new LinkedList<>();
                    energies.add(p.energyBreakdown.energyByResidue.get(p.sequence.indexOf(r)));
                    allEnergiesByAminoAcid.put(r.aminoAcid, energies);
                }
            }
        }

        Map<AminoAcid, Double> averageEnergyByAminoAcid = new HashMap<>();
        for (AminoAcid aminoAcid : allEnergiesByAminoAcid.keySet())
        {
            List<Double> energies = allEnergiesByAminoAcid.get(aminoAcid);
            double sum = 0.0;  // (could be issues with double overflowing)
            for (Double d : energies)
                sum += d;
            double average = sum / energies.size();
            averageEnergyByAminoAcid.put(aminoAcid, average);
        }

        return averageEnergyByAminoAcid;
    }

    /** 
    * A method to make the individual calls needed to run the process of finding reference energies
    * @return a map from amino acids to reference energies
    */
    public static Map<AminoAcid, Double> calculateReferenceEnergies()
    {
        // Draw random beta sheet backbones
        // Generate 5000 backbones and pick 500 lowest energy ones
        
        List<Peptide> sheets = BetaSheetGenerator.generateSheets(5, 1000, 5000, .01);
        Collections.sort(sheets);
        List<Peptide> lowEnergyBackbones = sheets;
        
        // Pick lowest energy structures
        for (int i = 0; i < 500; i++)
            lowEnergyBackbones.add(sheets.get(i));
        
        // Perform mutations
        List<Peptide> randomPeptides = generateRandomPeptides(lowEnergyBackbones);

        // Launch Monte Carlo jobs
        List<Peptide> minimizedRandomPeptides = new LinkedList<>();
        for (Peptide p : randomPeptides)
        {
            List<Peptide> monteCarloOutput = runMonteCarlo(p);
            // Minimize with AMOEBA
            List<Peptide> lowestEnergyMinimizedPeptides = minimize(monteCarloOutput);
            for (Peptide p2 : lowestEnergyMinimizedPeptides)
                minimizedRandomPeptides.add(p2);
        }

        // Breakdown energy for each peptide and average energies
        Map<AminoAcid, Double> referenceEnergies = averageEnergies(minimizedRandomPeptides);
        return referenceEnergies;
    }

    /** for testing */
    public static void main(String[] args)
    {
        // Artur, here is example code for how to rotamer pack:

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
        new GaussianInputFile(peptide).write("test_peptides/original.gjf");

        // create a list of pairs of backbone atoms to check for clashes
        // these pairs are defined as pairs of atoms that are more than 
        // this is kind of an expensive operation, but fortunately, it
        // only has to be done once per peptide composition, since the atom indices are not expected to change
        List<Pair<Integer,Integer>> backbonePairs = peptide.getBackbonePairs();

        // check peptide for self clashes
        if ( peptide.hasBackboneClash(backbonePairs) )
            {
                System.out.println("backbone clashes -- quitting");
                System.exit(1);
            }

        // create the A* iterator
        //
        // we must call this in a try-catch clause because it's possible the peptide could be in a bad conformation
        // in which it is impossible to place any rotamers
        FixedSequenceRotamerSpace fixedSequenceRotamerSpace = null;
        try
            {
                // note that the includeHN should be set to false if we want to do A* here
                fixedSequenceRotamerSpace = new FixedSequenceRotamerSpace(peptide, false);
            }
        catch (Exception e)
            {
                e.printStackTrace();
                System.exit(1);
            }

        // perform A* iteration, checking to see if the predicted and actual energies are the same
        AStarEnergyCalculator calculator = AStarEnergyCalculator.analyze(fixedSequenceRotamerSpace);
        RotamerIterator iterator = new RotamerIterator(fixedSequenceRotamerSpace.rotamerSpace, calculator.rotamerSelfEnergies, calculator.rotamerInteractionEnergies, 1000);
        List<RotamerIterator.Node> solutions = iterator.iterate();
        for (RotamerIterator.Node node : solutions)
            {
                List<Rotamer> rotamers = node.rotamers;
                Peptide thisPeptide = Rotamer.reconstitute(peptide, rotamers);
                // you can then minimize this peptide
                break;
                // you can iterate up to a maximum number of these solutions for minimization
            }
    }
}
