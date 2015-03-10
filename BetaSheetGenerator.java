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
 * This class takes a peptide of the form [Gly]n - D-Pro - Gly - [Gly]n and creates beta sheet conformations.
 *
 * The procedure is:
 *
 * 0. Create a peptide of the desired length.
 *
 * 1. Create a SheetUnit.  Every iteration, draw a random set of backbone angles for
 *    one of the residues.  Use X-ray-like angles for the hairpin and sheet-like (-120,+120) angles for the arms.
 *
 * 2. Check for very close clashes.  If they are not present, minimize on the OPLS forcefield.
 *
 * 3. Count up how many hydrogen bonds are present.  If this is an n-mer, we expect floor(n/2) hydrogen bonds.  If we
 *    don't see at least floor(n/2)-1 hydrogen bonds, throw out the structure.
 *
 * 4. Repeat steps 1-4 a certain number of times, discarding duplicates.  Return the lowest energy results.
 *
 * These steps run in parallel.  In practice, duplicate removal is accomplished by using a Map where the keys are
 * BackboneFingerprints.  These are basically lists of the backbone angles, where the angles are rounded to a certain
 * granularity.
 */
public class BetaSheetGenerator
{
    /** not instantiable */
    private BetaSheetGenerator()
    {
        throw new IllegalArgumentException("cannot be instantiated!");
    }

    /**
     * Creates beta sheets using all local cores.
     * @param armLength how many residues should go in each arm (this does not include the hairpin)
     * @param maxIterations how many Monte Carlo iterations to do per simulation
     * @param desiredNumberOfResults if this number of unique results is achieved (in aggregate), all runs will abort
     * @param deltaAlpha the amount to decrease the temperature by on every Monte Carlo iteration
     * @return the lowest energy beta sheet conformations
     */
    public static List<Peptide> generateSheets(int armLength, int maxIterations, int desiredNumberOfResults, double deltaAlpha)
    {
        // check the input
        if ( armLength < 2 )
            throw new IllegalArgumentException("must have at least two residues in each arm");
        if ( maxIterations < 1 )
            throw new IllegalArgumentException("you must do at least one iteration");
        if ( desiredNumberOfResults < 1 )
            throw new IllegalArgumentException("you must want more than one result");

        // create the seed peptide
        List<String> arm = new ArrayList<>(armLength);
        for (int i=0; i < armLength; i++)
            arm.add("standard_glycine");
        List<String> peptideStringList = new ArrayList<>(armLength*2 + 2);
        for (String s : arm)
            peptideStringList.add(s);
        peptideStringList.add("d_proline");
        peptideStringList.add("gly");
        for (String s : arm)
            peptideStringList.add(s);
        List<ProtoAminoAcid> sequence = ProtoAminoAcidDatabase.getSpecificSequence(peptideStringList);
        Peptide seedPeptide = PeptideFactory.createPeptide(sequence);
        seedPeptide = PeptideFactory.setHairpinAngles(seedPeptide);
        
        // queue the calculations
        Map<BackboneFingerprint,Peptide> resultMap = new ConcurrentHashMap<>();
        List<Future<Result>> futures = new ArrayList<>();
        for (int i=0; i < Settings.NUMBER_OF_THREADS; i++)
        //for (int i=0; i < 1; i++)
            {
                SheetUnit job = new SheetUnit(seedPeptide, resultMap, armLength, maxIterations, desiredNumberOfResults, deltaAlpha);
                Future<Result> f = GeneralThreadService.submit(job);
                futures.add(f);
            }

        // wait for them to finish
        GeneralThreadService.silentWaitForFutures(futures);

        // return the results
        return new ArrayList<Peptide>(resultMap.values());
    }

    /**
     * This runs a Monte Carlo simulation that tries to find low energy beta sheet conformations of a hairpin.
     */
    private static class SheetUnit implements WorkUnit
    {
        public final Peptide seedPeptide;
        public final Map<BackboneFingerprint,Peptide> resultMap;
        public final int armLength;
        public final int maxIterations;
        public final int desiredNumberOfResults;
        public final double deltaAlpha;
        public final int ID;

        public static final AtomicInteger ID_GENERATOR = new AtomicInteger();

        public SheetUnit(Peptide seedPeptide, Map<BackboneFingerprint,Peptide> resultMap, int armLength,
                         int maxIterations, int desiredNumberOfResults, double deltaAlpha)
        {
            this.seedPeptide = seedPeptide;
            this.resultMap = resultMap;
            this.armLength = armLength;
            this.maxIterations = maxIterations;
            this.desiredNumberOfResults = desiredNumberOfResults;
            this.deltaAlpha = deltaAlpha;
            ID = ID_GENERATOR.incrementAndGet();
        }

        public Result call()
        {
            // get some information before we start
            int numberOfResidues = seedPeptide.sequence.size();
            List<Integer> hairpinIndices = ImmutableList.of(numberOfResidues/2 - 1, numberOfResidues/2);
            int prolineResidueIndex = numberOfResidues/2 - 1;
            int glycineResidueIndex = prolineResidueIndex + 1;
            Peptide currentPeptide = seedPeptide;
            Peptide lastPeptide = seedPeptide;
            NormalDistribution prolinePhiDistribution   = new NormalDistribution(53.0, 5.0);
            NormalDistribution prolinePsiDistribution   = new NormalDistribution(-132.0, 5.0);
            NormalDistribution glycinePhiDistribution   = new NormalDistribution(-96.0, 15.0);
            NormalDistribution glycinePsiDistribution   = new NormalDistribution(9.0, 15.0);
            NormalDistribution betaNegativeDistribution = new NormalDistribution(-120.0, 15.0); 
            NormalDistribution betaPositiveDistribution = new NormalDistribution(120.0, 15.0);
            NormalDistribution omegaDistribution        = new NormalDistribution(180.0, 5.0);

            // do Monte Carlo
            double currentAlpha = 0.0;
            for (int iteration=1; iteration <= maxIterations; iteration++)
            {
                // quit if we have obtained the desired number of results
                if ( resultMap.size() >= desiredNumberOfResults )
                    {
                        System.out.printf("[%2d] Reached the maximum number of results.\n", ID);
                        break;
                    }

                // pick a residue at random
                int randomIndex = ThreadLocalRandom.current().nextInt(numberOfResidues);

                // make the mutation
                Peptide candidatePeptide = currentPeptide;
                if ( hairpinIndices.contains(randomIndex) )
                    {
                        // this is a hairpin position
                        candidatePeptide = BackboneMutator.setOmega (candidatePeptide, prolineResidueIndex, angleModulus(omegaDistribution.sample()));
                        candidatePeptide = BackboneMutator.setPhiPsi(candidatePeptide, prolineResidueIndex,
                                                                     angleModulus(prolinePhiDistribution.sample()),
                                                                     angleModulus(prolinePsiDistribution.sample()));
                        candidatePeptide = RotamerMutator.setChis   (candidatePeptide, prolineResidueIndex, ImmutableList.of(-22.0, 32.0, -29.0));
                        candidatePeptide = BackboneMutator.setOmega (candidatePeptide, glycineResidueIndex, angleModulus(omegaDistribution.sample()));
                        candidatePeptide = BackboneMutator.setPhiPsi(candidatePeptide, glycineResidueIndex,
                                                                     angleModulus(glycinePhiDistribution.sample()),
                                                                     angleModulus(glycinePsiDistribution.sample()));
                    }
                else
                    {
                        // this is a normal position
                        candidatePeptide = BackboneMutator.setOmega (candidatePeptide, randomIndex, angleModulus(omegaDistribution.sample()));
                        candidatePeptide = BackboneMutator.setPhiPsi(candidatePeptide, randomIndex,
                                                                     angleModulus(betaNegativeDistribution.sample()),
                                                                     angleModulus(betaPositiveDistribution.sample()));
                    }
                    
                // OPLS minimization
                //candidatePeptide = TinkerJob.minimize(candidatePeptide, Forcefield.OPLS, 100, false, true, false, false, false);
                try
                    {
                        candidatePeptide = TinkerJob.minimize(candidatePeptide, Forcefield.AMOEBA, 250, false, true, false, false, false);
                    }
                catch (Exception e)
                    {
                        System.out.printf("[%2d] %d of %d: rejected ( tinker error )\n", ID, iteration, maxIterations);
                        continue;
                    }

                // check for hydrogen bonds
                boolean isSheet = isSheet(candidatePeptide,1);
                if ( isSheet )
                    {
                        // only add to result map if this is a sheet
                        BackboneFingerprint backboneFingerprint = new BackboneFingerprint(candidatePeptide);
                        resultMap.put(backboneFingerprint, candidatePeptide);
                    }
                else
                    {
                        System.out.printf("[%2d] %d of %d: auto-rejected ( candidateE = %.2f, currentE = %.2f, isSheet = %b )\n", ID,
                                      iteration, maxIterations, candidatePeptide.energyBreakdown.totalEnergy, currentPeptide.energyBreakdown.totalEnergy, isSheet );
                        continue;
                    }
 
                // accept or reject
                boolean accept = MonteCarloJob.acceptChange(currentPeptide, candidatePeptide, currentAlpha);
                if ( accept )
                    {
                        currentPeptide = candidatePeptide;
                        System.out.printf("[%2d] %d of %d: accepted ( candidateE = %.2f, currentE = %.2f, isSheet = %b )\n", ID,
                                      iteration, maxIterations, candidatePeptide.energyBreakdown.totalEnergy, currentPeptide.energyBreakdown.totalEnergy, isSheet );
                    }
                else
                    System.out.printf("[%2d] %d of %d: rejected ( candidateE = %.2f, currentE = %.2f, isSheet = %b )\n", ID,
                                      iteration, maxIterations, candidatePeptide.energyBreakdown.totalEnergy, currentPeptide.energyBreakdown.totalEnergy, isSheet );
               
              }
            
            System.out.printf("[%2d] Reached the maximum number of iterations.\n", ID);
            return null;
        }
    }

    /**
     * Takes an angle and restricts it to the range [-180.0, 180.0] using the modulus.
     * @param d an angle
     * @return the restricted angle
     */
    public static double angleModulus(double d) 
    {
        double m = d % 360.0d;
        if (m < -180.0) m = m + 360.0;
        if (m >  180.0) m = m - 360.0;
        return m;
    }

    /**
     * Determines if the specified peptide is in a beta sheet conformation.
     * @param candidatePeptide the peptide to check
     * @param offset if the ideal number of interstrand H-bonds is n, then we ensure that there are at least n-offset bonds
     * @return true if this is in a beta sheet conformation
     */
    public static boolean isSheet(Peptide candidatePeptide, int offset)
    {
        // check for hydrogen bonds
        int numberOfResidues = candidatePeptide.sequence.size();
        int numberOfHydrogenBonds = 0;
        for (int i=0; i < (numberOfResidues / 2)-1; i++)
            {
                int j = numberOfResidues-i-1;
                Residue residueI = candidatePeptide.sequence.get(i);
                Residue residueJ = candidatePeptide.sequence.get(j);
                if ( residueJ.HN != null &&
                     Molecule.getDistance(residueI.O, residueJ.HN) < 2.5 && Molecule.getAngle(residueI.O, residueJ.HN, residueJ.N) > 120.0 )
                    numberOfHydrogenBonds++;
                if ( residueI.HN != null &&
                     Molecule.getDistance(residueJ.O, residueI.HN) < 2.5 && Molecule.getAngle(residueJ.O, residueI.HN, residueI.N) > 120.0 )
                    numberOfHydrogenBonds++;
            }
        if ( numberOfHydrogenBonds < (numberOfResidues / 2) - 1 - offset )
            return false;
        return true;
    }
    
    /**
     * Minimizes the specified sheets and checks them for whether they are beta sheets.
     * Does not use any solvation at all.
     * @param sheets the peptides to minimize
     * @param iterations the maximum number of iterations to use
     * @param forcefield which forcefield to use
     * @return the minimized beta sheetes
     */
    public static List<Peptide> minimizeSheets(List<Peptide> sheets, int iterations, Forcefield forcefield)
    {
        // run minimizations
        List<Future<Result>> futures = new ArrayList<>(sheets.size());
        for (Peptide p : sheets)
            {
                TinkerJob job = new TinkerJob(p, forcefield, iterations, false, false, false, false, false);
                Future<Result> f = GeneralThreadService.submit(job);
                futures.add(f);
            }
        //GeneralThreadService.silentWaitForFutures(futures);
        GeneralThreadService.waitForFutures(futures);

        // retrieve results
        List<Peptide> returnList = new ArrayList<>();
        for (Future<Result> f : futures)
            {
                try
                    {
                        TinkerJob.TinkerResult result = (TinkerJob.TinkerResult)f.get();
                        Peptide resultPeptide = result.minimizedPeptide;
                        //if ( isSheet(resultPeptide,0) )
                            returnList.add(resultPeptide);
                    }
                catch (Exception e) {}
            }
        return ImmutableList.copyOf(returnList);
    }
    
    /**
     * Minimizes sheets in serial.
     */
    public static List<Peptide> minimizeSheetsInSerial(List<Peptide> sheets, int iterations, Forcefield forcefield)
    {
        List<Peptide> results = new ArrayList<>();
        for (Peptide p : sheets)
            {
                TinkerJob job = new TinkerJob(p, forcefield, iterations, false, false, false, false, false);
                try
                    {
                        Peptide resultPeptide = job.call().minimizedPeptide;
                        if ( isSheet(resultPeptide, 0) )
                            results.add(resultPeptide);
                    }
                catch (Exception e) {}
            }
        return results;
    }

    /** Creates some beta sheets. */
    public static void main(String[] args)
    {
        DatabaseLoader.go();
        List<Peptide> sheets = generateSheets(6, 10, 10000, 0.01);
        sheets = new ArrayList<>(sheets);
        Collections.sort(sheets);

        // remove duplicates
        int maxResults = 100;
        double RMSDthreshold = 2.0;
        List<Peptide> results = new ArrayList<>();
        results.add(sheets.get(0));
        for (int i=1; i < sheets.size(); i++)
            {
                System.out.printf("superimposing %d of %d\r", i+1, sheets.size());
                Peptide candidatePeptide = sheets.get(i);
                Superposition superposition = Superposition.superimpose(candidatePeptide, results);
                boolean accept = true;
                System.out.println(superposition.RMSDs);
                for (Double RMSD : superposition.RMSDs)
                    {
                        // reject this candidate if it's too similar to an existing peptide
                        if ( RMSD <= RMSDthreshold )
                            {
                                accept = false;
                                break;
                            }
                    }
                if ( accept )
                    results.add(candidatePeptide);
                if ( results.size() >= maxResults )
                    break;
            }
        System.out.println("\nfinished superimposing");

        for ( int i=0; i < results.size(); i++ )
            {
                Peptide p = sheets.get(i);
                GaussianInputFile f = new GaussianInputFile(p);
                String filename = String.format("test_peptides/sheet_%02d.gjf", i);
                f.write(filename);
                System.out.printf("Sheet %02d: E = %7.2f\n", i, p.energyBreakdown.totalEnergy);
                
                filename = String.format("test_peptides/sheet_%02d.chk", i);
                p.checkpoint(filename);
            }
    }
}
