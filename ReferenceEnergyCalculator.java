import java.util.*;
import java.util.concurrent.*;
import java.util.concurrent.atomic.*;
import com.google.common.collect.*;
import java.io.*;

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
    * This method will perform mutaitons to a list of beta sheet peptides to randomly generate a list of peptides. 
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
                int index = ProtoAminoAcidLibrary.KEYS.indexOf(a);
                List<ProtoAminoAcid> VALUES = ProtoAminoAcidDatabase.VALUES.get(index);
                for (ProtoAminoAcid paa : VALUES)
                    {
                        // reject transition states
                        // decide upon this
                        if (paa.r.description.toLowerCase().indexOf("hairpin") > -1 )
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
                targetPAA = mutationOutcomes.get(0);
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
    public static List<Peptide> runMonteCarlo(Peptide randomPetptide)
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
            TinkerJob job = new TinkerJob(peptide, Forcefield.AMOEBA, 2000, false, true, true, false, true);
            TinkerResult result = job.call();
            Peptide newPeptide = result.minimizedPeptide; 
            
            PepideFingerprint fingerprint = new PeptideFingerprint(newPeptide, newPeptide.energyBreakdown.totalEnergy);
            
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
                    energies.add(p.energyBreakdown.energyByResidue.get(p.sequence.indexOf(r));
                    allEnergiesByAminoAcid.put(r.aminoAcid, energies);
                }
                else
                {
                    List<Double> energies = new LinkedList<>();
                    energies.add(p.energyBreakdown.energyByResidue.get(p.sequence.indexOf(r));
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
            for (Peptide p : lowestEnergyMinimizedPeptides)
                minimizedRandomPeptides.add(p);
        }

        // Breakdown energy for each peptide and average energies
        Map<AminoAcid, Double> referenceEnergies = averageEnergies(minimizedRandomPeptides);
    }

}
