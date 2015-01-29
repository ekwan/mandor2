import java.util.*;
import java.util.concurrent.*;
import java.util.concurrent.atomic.*;
import com.google.common.collect.*;
import org.apache.commons.io.*;
import java.io.*;

/**
 * This class collects together some methods to minimize and analyzepeptides.
 */
public class PeptideUtilities
{
    private PeptideUtilities()
    {
        throw new IllegalArgumentException("not instantiable");
    }

    /**
     * Minimizes a peptide on the OPLS forcefield in the gas phase.  An analysis is performed.
     * If the job fails, a null will be returned.
     * @param peptide the peptide to minimize
     * @param maxIterations the maximum number of iterations to perform during the minimization
     * @return the minimized peptide
     */
    public static Peptide minimizeWithOPLS(Peptide peptide, int maxIterations)
    {
        // check for sensible input
        if ( peptide == null )
            throw new NullPointerException("null input peptide is not allowed");
        if ( maxIterations < 1 )
            throw new IllegalArgumentException("positive number expected for maximum number of iterations");
        
        // perform minimization
        TinkerMinimizationJob job = new TinkerMinimizationJob(peptide, Forcefield.OPLS, "maxiter " + maxIterations + "\n\n");
        TinkerMinimizationJob.TinkerMinimizationResult result = job.call();
        Molecule minimizedMolecule = result.tinkerXYZOutputFile.molecule;
        Peptide newPeptide = peptide.setPositions(minimizedMolecule);

        System.out.println("grad is " + result.tinkerMinimizationLogFile.gradient);
        System.out.println(result.tinkerMinimizationLogFile.iterations + " iterations performed");
        
        // perform analysis
        TinkerAnalysisJob job2 = new TinkerAnalysisJob(newPeptide, Forcefield.OPLS, "\n\n");
        TinkerAnalysisJob.TinkerAnalysisResult result2 = job2.call();
        TinkerAnalyzeOutputFile outputFile2 = result2.tinkerAnalysisFile;
        
        System.out.println("OPLS energy by residue: " + outputFile2.energyByResidue);
        System.out.println("OPLS total energy is:   " + outputFile2.totalEnergy);

        // return new peptide
        EnergyBreakdown energyBreakdown = new EnergyBreakdown(outputFile2.energyByResidue, outputFile2.totalEnergy,
                                                              0.0, outputFile2.totalEnergy, null, EnergyBreakdown.Type.OPLS);
        newPeptide = newPeptide.setEnergyBreakdown(energyBreakdown);
        return newPeptide;
    }

    /**
     * Minimizes a peptide on the AMOEBA forcefield.  An analysis with solvation is performed.
     * If the job fails, a null will be returned.
     * @param peptide the peptide to minimize
     * @return the minimized peptide
     */
    public static Peptide minimizeWithAMOEBA(Peptide peptide)
    {
        return null;
    }

    /**
     * Minimizes multiple peptides on the OPLS forcefield.  Analyses are performed.
     * If a job fails, it will not appear in the final list.
     * @param peptides the peptides to minimize
     * @return the minimized peptides
     */
    public static List<Peptide> minimizeMultipleWithOPLS(List<Peptide> peptides)
    {
        return null;
    }

    /**
     * Minimizes multiple peptides on the OPLS forcefield.  Analyses are performed.
     * If a job fails, it will not appear in the final list.
     * @param peptides the peptides to minimize
     * @return the minimized peptides
     */
    public static List<Peptide> minimizeMultipleWithAMOEBA(List<Peptide> peptides)
    {
        return null;
    }

    /**
     * Takes a peptide with an energy breakdown and adds solvation to it.
     * @param peptide the peptide to analyze
     * @return the same peptide with an EnergyBreakdown that accounts for solvation
     */
    public static Peptide analyzeWithSolvation(Peptide peptide) 
    {
        SolvationUnit unit = new SolvationUnit(peptide);
        return unit.call().peptide;
    }

    /**
     * Takes a list of peptides with energy breakdowns and adds solvation energies.
     * @param peptides the peptides to analyze
     * @return the same peptides with EnergyBreakdowns that account for solvation
     */
    public static List<Peptide> analyzeWithSolvation(List<Peptide> peptides)
    {
        List<Future<Result>> futures = new ArrayList<>(peptides.size());
        for (Peptide p : peptides)
            {
                SolvationUnit unit = new SolvationUnit(p);
                Future<Result> f = GeneralThreadService.submit(unit);
                futures.add(f);
            }
        GeneralThreadService.silentWaitForFutures(futures);
        List<Peptide> results = new ArrayList<>(peptides.size());
        for (Future<Result> f : futures)
            {
                try
                    {
                        SolvationResult result = (SolvationResult)f.get();
                        results.add(result.peptide);
                    }
                catch (Exception e)
                    {
                        e.printStackTrace();
                    }
            }
        return ImmutableList.copyOf(results);
    }

    public static class SolvationUnit implements WorkUnit
    {
        /** The peptide to calculate a solvation energy for. */
        public final Peptide peptide;

        /** Creates a SolvationUnit. */
        public SolvationUnit(Peptide peptide)
        {
            if ( peptide.energyBreakdown == null || peptide.energyBreakdown == EnergyBreakdown.BLANK )
                throw new IllegalArgumentException("expected filled out energy breakdown first");
            if ( peptide.energyBreakdown.solvationEnergy != 0.0 )
                throw new IllegalArgumentException("solvation energy already set");
            this.peptide = peptide;
        }

        @Override
        public SolvationResult call()
        {
            // compute solvation
            List<Double> SASAlist = DCLMAreaCalculator.INSTANCE.calculateSASA(peptide);
            List<Double> energies = new ArrayList<>(SASAlist.size());
            double solvationEnergy = 0.0;
            for (int i=0; i < SASAlist.size(); i++)
                {
                    double surfaceArea = SASAlist.get(i);
                    double surfaceTension = peptide.contents.get(i).surfaceTension;
                    double energy = surfaceArea * surfaceTension;
                    energies.add(energy);
                    //System.out.printf("%3d  %8.2f  %8.2f\n", i+1, surfaceArea, energy);
                    solvationEnergy += energy;
                }
            //System.out.printf("solvation energy: %8.2f\n", solvationEnergy);

            // write out debug file for OMNISOL
            //String omnisolString = "SM5.0R\n& IOFR=1.4459 ALPHA=0.15 BETA=0.02 GAMMA=38.39\n& FACARB=0.00 FEHALO=0.75 SOLVNT=GENORG\npeptide (solvent : chloroform)\n\n";
            //omnisolString += newPeptide.toOmnisolString();
            //InputFileFormat.writeStringToDisk(omnisolString, "debug.dat");

            // figure out which energies belong to which residues
            double[] solvationEnergiesByResidue = new double[peptide.sequence.size()]; 
            int numberOfAtoms = 0; 
            for (int i=0; i < peptide.sequence.size(); i++)
                {
                    Residue residue = peptide.sequence.get(i);
                    for (Atom a : residue.atoms)
                        {
                            numberOfAtoms++;
                            int atomIndex = peptide.contents.indexOf(a);
                            if ( atomIndex == -1 )
                                throw new IllegalArgumentException("atom not in molecule");
                            solvationEnergiesByResidue[i] += energies.get(atomIndex);
                        }
                }
            System.out.println(Arrays.toString(solvationEnergiesByResidue));
            if (numberOfAtoms != peptide.contents.size())
                throw new IllegalArgumentException("failed to process solvation energies for all atoms");

            // update energy breakdown
            EnergyBreakdown oldEnergyBreakdown = peptide.energyBreakdown;

            List<Double> energyByResidue = new ArrayList<>(oldEnergyBreakdown.energyByResidue);
            System.out.println("old: " + energyByResidue.toString());
            System.out.println("+:   " + Arrays.toString(solvationEnergiesByResidue));
            for (int i=0; i < energyByResidue.size(); i++)
                {
                    double oldEnergy = energyByResidue.get(i);
                    double newEnergy = oldEnergy + solvationEnergiesByResidue[i];
                    energyByResidue.set(i, newEnergy);
                }
            energyByResidue = ImmutableList.copyOf(energyByResidue);
            System.out.println("new: " + energyByResidue.toString());

            double totalEnergy = oldEnergyBreakdown.totalEnergy + solvationEnergy;
            double potentialEnergy = oldEnergyBreakdown.potentialEnergy;
            
            EnergyBreakdown newEnergyBreakdown = new EnergyBreakdown(energyByResidue, totalEnergy, solvationEnergy,
                                                                     potentialEnergy, SASAlist, oldEnergyBreakdown.type);

            // return result
            Peptide newPeptide = peptide.setEnergyBreakdown(newEnergyBreakdown);
            return new SolvationResult(newPeptide);
        }
    }

    /** This is a wrapper object that contains a peptide from a SolvationJob.  This peptide has solvation included in its EnergyBreakdown. */
    public static class SolvationResult implements Result
    {
        public final Peptide peptide;

        public SolvationResult(Peptide peptide)
        {
            this.peptide = peptide;
        }
    }

    /** For testing. */
    public static void main(String[] args)
    {
        DatabaseLoader.go();
        System.out.println("building");
        List<ProtoAminoAcid> sequence = ProtoAminoAcidDatabase.getSpecificSequence("arg","met","standard_ala","gly","d_proline", "gly", "phe", "val", "hd",         "l_pro");
        Peptide peptide = PeptideFactory.createPeptide(sequence);

        for (int i=0; i < peptide.sequence.size(); i++)
            {
                peptide = BackboneMutator.mutateOmega(peptide, i);
                peptide = BackboneMutator.mutatePhiPsi(peptide, i);
                peptide = RotamerMutator.mutateChis(peptide, i);
            }
        peptide = PeptideFactory.setHairpinAngles(peptide);
        System.out.println("done");

        Peptide newPeptide = minimizeWithOPLS(peptide, 1000);
        newPeptide = analyzeWithSolvation(newPeptide);
    }
}
