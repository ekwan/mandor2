import java.util.*;
import java.util.concurrent.*;
import java.util.concurrent.atomic.*;
import com.google.common.collect.*;
import java.io.*;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.geometry.euclidean.threed.*;
import org.jgrapht.*;
import org.jgrapht.graph.*;
import org.jgrapht.alg.*;

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
    public static final NUMBER_BETA_SHEETS = 500;

    /** The number of lowest energy structures to keep for each starting random peptide following fixed sequence Monte Carlo */
    public static final NUMBER_LOWEST_ENERGY_STRUCTURES_TO_KEEP = 10;

    /** The number of structures to minimize using AMOEBA following the fixed sequence Monte Carlo */
    public static final NUMBER_STRUCTURES_TO_MINIMIZE = 1000;

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
    public static List<Peptide> performMutations(List<Peptide> betaSheets)
    {

    }

    
    /** 
    * This method will run the fixed sequence Monte Carlo process.
    * @param randomPeptide a randomly generated peptide from performMutations
    * @return a list of size NUMBER_STRUCTURES_TO_MINIMIZE that will contain the lowest energy peptides from a fixed sequence simulation
    */
    public static List<Peptide> runMonteCarlo(Peptide randomPetptide)
    {

    }

    /**
    * A method for performing AMOEBA minimizations that also removes duplicates and returns the lowest energy structures only
    * @param monteCarloPeptides the peptides generated from the Monte Carlo process which includes NUMBER_STRUCTURES_TO_MINIMIZE peptides
    * @return the NUMBER_LOWEST_ENERGY_STRUCTURES_TO_KEEP peptides that will be used to find the average amino acid reference energy 
    */
    public static List<Peptide> minimize(List<Peptide> monteCarloPeptides)
    {

    }

    /** 
    * This method will calculate the average energy of each amino acid by taking an energy breakdown of each of the provided structures
    * @param structures a random set of low-energy structures 
    * @return a map from amino acids to average energies
    */
    public static Map<AminoAcid, Double> averageEnergies(List<Peptide> structures)
    {

    }

    /** 
    * A method to make the individual calls needed to run the process of finding reference energies
    * @return a map from amino acids to reference energies
    */
    public static Map<AminoAcid, Double> calculateReferenceEnergies()
    {

    }

}
