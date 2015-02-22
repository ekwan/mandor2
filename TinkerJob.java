import java.util.*;
import java.util.concurrent.*;
import java.util.concurrent.atomic.*;
import com.google.common.collect.*;
import org.apache.commons.io.*;
import java.io.*;

/**
 * This is an all purpose work unit for doing Tinker calculations with peptides.
 */
public class TinkerJob implements WorkUnit
{
    /** The peptide to minimize. */
    public final Peptide peptide;

    /** The forcefield to use. */
    public final Forcefield forcefield;

    /** The maximum number of minimization iterations to do.  Note that TinkerMinimizationJob.MAX_JOB_TIME also applies. */
    public final int maxIterations;

    /** Whether to use Tinker solvation during the minimization. */
    public final boolean solvateDuringMinimization;

    /** Whether we should add an approximate solvation correction at the end of the minimization.  Cannot be turned of if analysis is turned on. */
    public final boolean approximateSolvationSinglePoint;

    /** Whether to break down the energies into residue-by-residue terms. */
    public final boolean doAnalysis;

    /** Whether to use Tinker solvation during the analysis, if any.  If there's no analysis, this is ignored. */
    public final boolean solvateDuringAnalyze;

    /** Whether to do solvation with approximate OMNISOL surface tensions.  If there's no analysis, this is ignored. */
    public final boolean approximateSolvationDuringAnalyze;

    /** See field comments for restrictions. */
    public TinkerJob(Peptide peptide, Forcefield forcefield, int maxIterations, boolean solvateDuringMinimization,
                     boolean approximateSolvationSinglePoint, boolean doAnalysis, boolean solvateDuringAnalyze,
                     boolean approximateSolvationDuringAnalyze)
    {
        if ( peptide == null )
            throw new NullPointerException("null peptide not allowed");
        this.peptide = peptide;

        if ( forcefield == null )
            throw new NullPointerException("must specify a forcefield");
        this.forcefield = forcefield;

        if ( maxIterations < 1 )
            throw new IllegalArgumentException("must have at least one iteration");
        this.maxIterations = maxIterations;

        if ( solvateDuringMinimization && approximateSolvationSinglePoint)
            throw new IllegalArgumentException("approximate solvation single point is redundant");
        if ( doAnalysis && approximateSolvationSinglePoint )
            throw new IllegalArgumentException("can't do approximate solvation single point without a full analysis");
        if ( doAnalysis && solvateDuringAnalyze && approximateSolvationDuringAnalyze)
            throw new IllegalArgumentException("can't do two kinds of solvation at once");
        if ( doAnalysis && forcefield == Forcefield.AMOEBA && solvateDuringAnalyze )
            throw new IllegalArgumentException("unfortunately an analysis of AMOEBA-GK solvation does not work yet");

        this.solvateDuringMinimization = solvateDuringMinimization;
        this.approximateSolvationSinglePoint = approximateSolvationSinglePoint;
        this.doAnalysis = doAnalysis;
        this.solvateDuringAnalyze = solvateDuringAnalyze;
        this.approximateSolvationDuringAnalyze = approximateSolvationDuringAnalyze;
    }

    /** Performs the composite tinker jobs. */
    public TinkerResult call()
    {
        // set up parameters
        String solvationString = null;
        if ( forcefield == Forcefield.OPLS )
            solvationString = "solvate gb\n";
        else if ( forcefield == Forcefield.AMOEBA )
            solvationString = "solvate gk\n";
        else
            throw new IllegalArgumentException("unknown forcefield");
        solvationString += "\n";
        
        // perform minimization
        String parameters = String.format("maxiter %d\n", maxIterations);
        if ( solvateDuringMinimization )
            parameters += solvationString;
        parameters += "\n";
        TinkerMinimizationJob job = new TinkerMinimizationJob(peptide, forcefield, parameters);
        TinkerMinimizationJob.TinkerMinimizationResult result = job.call();
        Molecule minimizedMolecule = result.tinkerXYZOutputFile.molecule;
        Peptide newPeptide = peptide.setPositions(minimizedMolecule);

        // deal with energies
        TinkerMinimizationJob.TinkerMinimizationLogFile tinkerMinimizationLogFile = result.tinkerMinimizationLogFile;
        double potentialEnergy = tinkerMinimizationLogFile.energy;
        double gradient = tinkerMinimizationLogFile.gradient;
        int iterations = tinkerMinimizationLogFile.iterations;

        // return result now if nothing else was requested
        if ( !solvateDuringMinimization && !approximateSolvationSinglePoint && !doAnalysis )
            {
                EnergyBreakdown energyBreakdown = new EnergyBreakdown(null, potentialEnergy, 0.0, potentialEnergy, null, forcefield);
                newPeptide = newPeptide.setEnergyBreakdown(energyBreakdown);
                return new TinkerResult(newPeptide);
            }

        // debugging
        //System.out.printf("Forcefield:       %s\n", forcefield.toString());
        //System.out.printf("Iterations:       %d\n", iterations);
        //System.out.printf("RMS Gradient:     %.2f\n", gradient);
        //System.out.printf("Potential Energy: %.2f\n", potentialEnergy);

        // add approximate solvation single point if requested
        if ( approximateSolvationSinglePoint )
            {
                double solvationEnergy = 0.0;
                List<Double> SASAlist = null;
                try { SASAlist = new DCLMAreaCalculator(0.0).calculateSASA(newPeptide); }
                catch (Exception e) { e.printStackTrace(); SASAlist = ShrakeRupleyCalculator.INSTANCE.calculateSASA(newPeptide); }
                for (int i=0; i < SASAlist.size(); i++)
                    {
                        double surfaceArea = SASAlist.get(i);
                        double surfaceTension = peptide.contents.get(i).surfaceTension;
                        double energy = surfaceArea * surfaceTension;
                        //System.out.printf("%3d  %8.2f  %8.2f\n", i+1, surfaceArea, energy);
                        solvationEnergy += energy;
                    }
                double totalEnergy = potentialEnergy + solvationEnergy;
                //System.out.printf("Solvation Energy: %.2f\n", solvationEnergy);
                //System.out.printf("Total Energy:     %.2f\n", totalEnergy);
                EnergyBreakdown energyBreakdown = new EnergyBreakdown(null, totalEnergy, solvationEnergy, potentialEnergy, SASAlist, forcefield);
                return new TinkerResult(newPeptide.setEnergyBreakdown(energyBreakdown));
            }
       
        // perform analysis if requested
        // note that if tinker solvation is used, the solvation energy will show up as 0.0 in the energy breakdown because it's assumed
        // to be part of the potential energy
        parameters = "\n";
        if ( solvateDuringMinimization )
            parameters += solvationString;
        parameters += "\n";
        TinkerAnalysisJob job2 = new TinkerAnalysisJob(newPeptide, forcefield, "\n\n");
        TinkerAnalysisJob.TinkerAnalysisResult result2 = job2.call();
        TinkerAnalyzeOutputFile outputFile2 = result2.tinkerAnalysisFile;
        EnergyBreakdown energyBreakdown = new EnergyBreakdown(outputFile2.energyByResidue, outputFile2.totalEnergy,
                                                              0.0, outputFile2.totalEnergy, null, forcefield);
        newPeptide = newPeptide.setEnergyBreakdown(energyBreakdown);

        // return the peptide with an energy breakdown if we are finished
        if ( !approximateSolvationDuringAnalyze )
            return new TinkerResult(newPeptide);
        else
            {
                // otherwise, compute approximate solvation
                List<Double> SASAlist = null;
                try { SASAlist = new DCLMAreaCalculator(0.0).calculateSASA(newPeptide); }
                catch (Exception e) { e.printStackTrace(); SASAlist = ShrakeRupleyCalculator.INSTANCE.calculateSASA(newPeptide); }
                List<Double> energies = new ArrayList<>(SASAlist.size());
                double solvationEnergy = 0.0;
                for (int i=0; i < SASAlist.size(); i++)
                    {
                        double surfaceArea = SASAlist.get(i);
                        double surfaceTension = newPeptide.contents.get(i).surfaceTension;
                        double energy = surfaceArea * surfaceTension;
                        energies.add(energy);
                        //System.out.printf("%3d  %8.2f  %8.2f\n", i+1, surfaceArea, energy);
                        solvationEnergy += energy;
                    }
                //System.out.printf("Solvation energy: %.2f\n", solvationEnergy);

                // write out debug file for OMNISOL
                //String omnisolString = "SM5.0R\n& IOFR=1.4459 ALPHA=0.15 BETA=0.02 GAMMA=38.39\n& FACARB=0.00 FEHALO=0.75 SOLVNT=GENORG\npeptide (solvent : chloroform)\n\n";
                //omnisolString += newPeptide.toOmnisolString();
                //InputFileFormat.writeStringToDisk(omnisolString, "debug.dat");

                // figure out which energies belong to which residues
                double[] solvationEnergiesByResidue = new double[newPeptide.sequence.size()]; 
                int numberOfAtoms = 0; 
                for (int i=0; i < newPeptide.sequence.size(); i++)
                    {
                        Residue residue = newPeptide.sequence.get(i);
                        for (Atom a : residue.atoms)
                            {
                                numberOfAtoms++;
                                int atomIndex = newPeptide.contents.indexOf(a);
                                if ( atomIndex == -1 )
                                    throw new IllegalArgumentException("atom not in molecule");
                                solvationEnergiesByResidue[i] += energies.get(atomIndex);
                            }
                    }
                if (numberOfAtoms != peptide.contents.size())
                    throw new IllegalArgumentException("failed to process solvation energies for all atoms");

                // update energy breakdown
                EnergyBreakdown oldEnergyBreakdown = newPeptide.energyBreakdown;

                List<Double> energyByResidue = new ArrayList<>(oldEnergyBreakdown.energyByResidue);
                //System.out.println("old: " + energyByResidue.toString());
                //System.out.println("+:   " + Arrays.toString(solvationEnergiesByResidue));
                for (int i=0; i < energyByResidue.size(); i++)
                    {
                        double oldEnergy = energyByResidue.get(i);
                        double newEnergy = oldEnergy + solvationEnergiesByResidue[i];
                        energyByResidue.set(i, newEnergy);
                    }
                energyByResidue = ImmutableList.copyOf(energyByResidue);
                //System.out.println("new: " + energyByResidue.toString());

                double totalEnergy = potentialEnergy + solvationEnergy;
                //System.out.printf("Solvation Energy: %.2f\n", solvationEnergy);
                //System.out.printf("Total Energy:     %.2f\n", totalEnergy);
                EnergyBreakdown newEnergyBreakdown = new EnergyBreakdown(energyByResidue, totalEnergy, solvationEnergy,
                                                                         potentialEnergy, SASAlist, forcefield);

                // return result
                return new TinkerResult(newPeptide.setEnergyBreakdown(newEnergyBreakdown));
            }
    }

    public static class TinkerResult implements Result
    {
        public final Peptide minimizedPeptide;

        public TinkerResult(Peptide minimizedPeptide)
        {
            this.minimizedPeptide = minimizedPeptide;
        }
    }

    /**
     * Minimizes a single peptide using the specified parameters.  Exceptions are not caught.
     */
    public static Peptide minimize(Peptide peptide, Forcefield forcefield, int maxIterations, boolean solvateDuringMinimization,
                                   boolean approximateSolvationSinglePoint, boolean doAnalysis, boolean solvateDuringAnalyze,
                                   boolean approximateSolvationDuringAnalyze)
    {
        TinkerJob job = new TinkerJob(peptide, forcefield, maxIterations, solvateDuringMinimization, approximateSolvationSinglePoint,
                                      doAnalysis, solvateDuringAnalyze, approximateSolvationDuringAnalyze);
        return job.call().minimizedPeptide;
    }

    /**
     * Minimizes multiple peptides using the specified parameters.  If a job runs into an exception, it will not be included in the final result.
     */
    public static List<Peptide> minimize(List<Peptide> peptides, Forcefield forcefield, int maxIterations, boolean solvateDuringMinimization,
                                         boolean approximateSolvationSinglePoint, boolean doAnalysis, boolean solvateDuringAnalyze,
                                         boolean approximateSolvationDuringAnalyze)
    {
        List<Future<Result>> futures = new ArrayList<>(peptides.size());
        for (Peptide p : peptides)
            {
                TinkerJob job = new TinkerJob(p, forcefield, maxIterations, solvateDuringMinimization,
                                              approximateSolvationSinglePoint, doAnalysis, solvateDuringAnalyze, approximateSolvationDuringAnalyze);
                Future<Result> f = GeneralThreadService.submit(job);
                futures.add(f);
            }
        GeneralThreadService.waitForFutures(futures);
        List<Peptide> results = new ArrayList<>(peptides.size());
        for (Future<Result> f : futures)
            {
                try
                    {
                        TinkerResult result = (TinkerResult)f.get();
                        results.add(result.minimizedPeptide);
                    }
                catch (Exception e)
                    {
                        //e.printStackTrace();
                    }
            }
        Collections.sort(results);
        return results;
    }

    /**
     * Performs an analysis and energy breakdown on all the specified peptides.
     * Analyses will be done in the gas phase
     * @param peptides the peptides to analyze
     * @param forcefield the forcefield to use
     * @return the analyzed peptides
     */
    public static List<Peptide> analyze(List<Peptide> peptides, Forcefield forcefield)
    {
        Map<Future<Result>,Peptide> futureMap = new HashMap<>();
        List<Future<Result>> futures = new ArrayList<>(peptides.size());
        for (Peptide p : peptides)
            {
                TinkerAnalysisJob job = new TinkerAnalysisJob(p, forcefield, "\n\n");
                Future<Result> f = GeneralThreadService.submit(job);
                futures.add(f);
                futureMap.put(f, p);
            }
        GeneralThreadService.silentWaitForFutures(futures);
        List<Peptide> results = new ArrayList<>(peptides.size());
        for (Future<Result> f : futures)
            {
                try
                    {
                        TinkerAnalysisJob.TinkerAnalysisResult result = (TinkerAnalysisJob.TinkerAnalysisResult)f.get();
                        TinkerAnalyzeOutputFile outputFile = result.tinkerAnalysisFile;
                        EnergyBreakdown energyBreakdown = new EnergyBreakdown(outputFile.energyByResidue, outputFile.totalEnergy,
                                                                              0.0, outputFile.totalEnergy, null, forcefield);
                        Peptide oldPeptide = futureMap.get(f);
                        Peptide newPeptide = oldPeptide.setEnergyBreakdown(energyBreakdown);
                        results.add(newPeptide);                    
                    }
                catch (Exception e) {}
            }
        return results;
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

        peptide = minimize(peptide, Forcefield.OPLS, 1000, // max iterations
                           false,  // do tinker solvation during the minimization
                           false,   // correct minimization energy with an approximate solvation energy
                           false,  // perform analysis
                           false,  // use tinker solvation during the analysis
                           true); // use approximate solvation during the analysis
    }
}
