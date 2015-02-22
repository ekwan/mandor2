import java.io.*;
import java.util.*;
import com.google.common.collect.*;
import java.util.concurrent.*;
import java.util.concurrent.atomic.*;
import org.apache.commons.math3.geometry.euclidean.threed.Vector3D;
import org.apache.commons.io.FileUtils;
import org.apache.commons.math3.distribution.NormalDistribution;

/**
 * This minimizes the conformation of one peptide in serial.
 */
public class FixedSequenceMonteCarloJob extends MonteCarloJob implements Serializable
{
    /** For serialization. */
    public static final long serialVersionUID = 1L;

    /** If a non-hairpin proline is selected for mutation, the chance of rotating it to the cis conformation. */
    public static final double CIS_PROLINE_PROBABILITY = 0.10;

    /** For setting omegas. */
    public static final NormalDistribution TRANS_DISTRIBUTION = new NormalDistribution(180.0,5.0);
    
    /** For setting omegas. */
    public static final NormalDistribution CIS_DISTRIBUTION = new NormalDistribution(0.0,5.0);
   
    /**
     * The allowed (phi,psi) angles for the hairpin D-proline and glycine.  Entries 0 and 1 are phi and psi for
     * the proline, respectively.  Entries 2 and 3 are phi and psi for the glycine, respectively.
     */
    public static final DiscreteProbabilityDistribution<List<Double>> HAIRPIN_DISTRIBUTION;

    /** Static initializer. */
    static
    { 
        List<List<Double>> outcomes = new ArrayList<>();
        List<Double> probabilities = new ArrayList<>();
        
        // type II' beta hairpin
        outcomes.add(ImmutableList.of(60.0,-120.0,-80.0,0.0));
        probabilities.add(0.70);

        // type I' beta hairpin
        outcomes.add(ImmutableList.of(60.0,30.0,90.0,0.0));
        probabilities.add(0.15);

        // gamma turn
        outcomes.add(ImmutableList.of(75.0,-64.0,-81.0,9.0));
        probabilities.add(0.15);

        HAIRPIN_DISTRIBUTION = new DiscreteProbabilityDistribution<List<Double>>(outcomes,probabilities);
    }
    
    /** A list of backbone atom pairs to check for clashes. */
    public final List<Pair<Integer, Integer>> backbonePairs;

    /** How many rotamers should be minimized on each iteration. */
    public final int rotamersPerIteration;

    /** Stores the best results. */
    public final MonteCarloJob.PeptideList bestPeptides;

    /** Where to save the result to periodically. */
    public final String checkpointFilename;

    /** Unique ID. */
    private long serverID;

    /**
     * Sets up job but doesn't run it.
     * @param startingPeptide the starting structure
     * @param deltaAlpha the cooling rate
     * @param maxIterations how many MC iterations to perform
     * @param maxSize the number of peptides to keep
     * @param how many rotamers should be minimized on each iteration
     * @param checkpointFilename
     */
    public FixedSequenceMonteCarloJob(Peptide startingPeptide, double deltaAlpha, int maxIterations, int maxSize, int rotamersPerIteration, String checkpointFilename)
    {
        super(startingPeptide, deltaAlpha, maxIterations);
        this.backbonePairs = startingPeptide.getBackbonePairs();
        this.bestPeptides = new MonteCarloJob.PeptideList(maxSize);
        this.rotamersPerIteration = rotamersPerIteration;
        this.checkpointFilename = checkpointFilename;
        this.serverID = RemoteWorkUnit.ID_GENERATOR.incrementAndGet();
    }

    /**
     * Randomly changes the backbone and then rotamer packs.  The procedure is:<p>
     * (1) Select a residue at random.<p>
     * (2) If this is a hairpin residue, set the hairpin to one of several pre-defined conformations.<p>
     * (3) If this is not a hairpin residue, choose (omega,phi,psi) for it.  If a proline was selected,
     *     give a cis-proline angle with some probability.<p>
     * (4) Rotamer pack and pick lowest energy candidate sturcture. <p>
     * Obviously, the sequence cannot be changed.
     * @param peptide the peptide whose backbone we should perturb
     * @return the perturbed peptide
     */
    public Peptide mutate(Peptide peptide)
    {
        // random number generator
        ThreadLocalRandom random = ThreadLocalRandom.current();

        // the new peptide
        Peptide newPeptide = peptide;
        FixedSequenceRotamerSpace fixedSequenceRotamerSpace = null;
        
        // choose a residue at random
        int randomIndex = -1;
        Residue residue = null;

        while (true)
            {
                randomIndex = random.nextInt(peptide.sequence.size());
                residue = peptide.sequence.get(randomIndex);
                
                // if this is a hairpin residue, give a high chance of rolling again
                if ( residue.isHairpin && random.nextDouble() > 0.10 )
                    continue;
                break;
            }

        // perturb the backbone until there's no backbone clash
        int microiteration = 0;
        while (true)
            {
                // reset
                newPeptide = peptide;
                microiteration++;

                // make the mutation
                if ( residue.isHairpin )
                    {
                        // if this is a hairpin residue, set the hairpin
                        // this is hard-coded for a 12-mer
                        
                        // choose the phi,psi angles
                        List<Double> hairpinAngles = HAIRPIN_DISTRIBUTION.getRandom();
                        double prolinePhiAngle = hairpinAngles.get(0);
                        double prolinePsiAngle = hairpinAngles.get(1);
                        double glycinePhiAngle = hairpinAngles.get(2);
                        double glycinePsiAngle = hairpinAngles.get(3);

                        System.out.printf("[%3d] Set hairpin backbone to: %.0f, %.0f, %.0f, %.0f on microiteration %d.\n", serverID, prolinePhiAngle, prolinePsiAngle, glycinePhiAngle, glycinePsiAngle, microiteration);

                        // choose the omega angles
                        double prolineOmegaAngle = BetaSheetGenerator.angleModulus(TRANS_DISTRIBUTION.sample());
                        double glycineOmegaAngle = BetaSheetGenerator.angleModulus(TRANS_DISTRIBUTION.sample());

                        // apply the changes
                        int forbiddenIndex = (newPeptide.sequence.size() / 2 ) - 1;
                        Residue proline = newPeptide.sequence.get(forbiddenIndex);
                        newPeptide = BackboneMutator.setPhiPsi(newPeptide, proline, prolinePhiAngle, prolinePsiAngle);
                        proline = newPeptide.sequence.get(forbiddenIndex);
                        newPeptide = BackboneMutator.setOmega(newPeptide, proline, prolineOmegaAngle);

                        Residue glycine = newPeptide.sequence.get(forbiddenIndex+1);
                        newPeptide = BackboneMutator.setPhiPsi(newPeptide, glycine, glycinePhiAngle, glycinePsiAngle);
                        glycine = newPeptide.sequence.get(forbiddenIndex+1);
                        newPeptide = BackboneMutator.setOmega(newPeptide, glycine, glycineOmegaAngle);
                    }
                else
                    {
                        // otherwise, choose a random (omega,phi,psi)
                        System.out.printf("[%3d] Setting backbone angles for residue %d (%s) on microiteration %d.\n", serverID, randomIndex, residue.aminoAcid.shortName, microiteration);
                        if ( residue.aminoAcid.isProline() )
                            {
                                // decide whether we should set this to cis
                                double roll = random.nextDouble();
                                double omegaAngle = 180.0;
                                if ( roll < CIS_PROLINE_PROBABILITY )
                                    omegaAngle = BetaSheetGenerator.angleModulus(CIS_DISTRIBUTION.sample()); 
                                else
                                    omegaAngle = BetaSheetGenerator.angleModulus(TRANS_DISTRIBUTION.sample());

                                // apply the change
                                newPeptide = BackboneMutator.setOmega(newPeptide, residue, omegaAngle);

                                // choose a random (phi,psi)
                                newPeptide = BackboneMutator.mutatePhiPsi(newPeptide,randomIndex+1);
                            }
                        else
                            {
                                // choose a random (omega,phi,psi)
                                newPeptide = BackboneMutator.mutateOmega(newPeptide,randomIndex);
                                newPeptide = BackboneMutator.mutatePhiPsi(newPeptide,randomIndex);
                            }
                    }
                    
                // check peptide for self clashes
                if ( peptide.hasBackboneClash(backbonePairs) )
                    {
                        System.out.printf("[%3d] rejected (backbone clash) on microiteration %d\n", serverID, microiteration);
                        continue;
                    }

                // rotamer pack with an A* iteration
                //
                // we must call this in a try-catch clause because it's possible the peptide could be in a bad conformation
                // in which it is impossible to place any rotamers
                try
                    {
                        // note that the includeHN should be set to false if we want to do A* here
                        fixedSequenceRotamerSpace = new FixedSequenceRotamerSpace(peptide, false);
                    }
                catch (Exception e)
                    {
                        if ( e instanceof IllegalArgumentException )
                            System.out.printf("[%3d] rejected (rotamer packing / failure %s) on microiteration %d\n", serverID, e.getMessage(), microiteration);
                        else
                            {
                                System.out.printf("[%3d] rejected (rotamer packing unknown error) on microiteration %d\n", serverID, microiteration);
                                e.printStackTrace();
                            }
                        continue;
                    }

                // perform A* iteration, checking to see if the predicted and actual energies are the same
                AStarEnergyCalculator calculator = AStarEnergyCalculator.analyze(fixedSequenceRotamerSpace);
                RotamerIterator iterator = new RotamerIterator(fixedSequenceRotamerSpace.rotamerSpace, calculator.rotamerSelfEnergies, calculator.rotamerInteractionEnergies, rotamersPerIteration);
                List<Peptide> results = new ArrayList<>(rotamersPerIteration+100);
                List<RotamerIterator.Node> solutions = iterator.iterate();
                if ( solutions.size() == 0 )
                    {
                        System.out.printf("[%3d] rejected (no A* solutions) on microiteration %d\n", serverID, microiteration);
                        continue;
                    }

                // minimize peptides in serial
                for (RotamerIterator.Node node : solutions)
                    {
                        // solutions come out lowest energy first
                        List<Rotamer> rotamers = node.rotamers;
                        Peptide thisPeptide = Rotamer.reconstitute(peptide, rotamers);

                        // minimize with OPLS with single point GB solvation, throw out peptides that error
                        Peptide minimizedPeptide = null;
                        try { minimizedPeptide = TinkerJob.minimize(thisPeptide, Forcefield.OPLS, 2000, false, false, true, true, false); }
                        catch (Exception e) { continue; }
                        results.add(minimizedPeptide);
                        bestPeptides.add(minimizedPeptide);
                        if ( results.size() >= rotamersPerIteration )
                            break;
                    }

                // check that we have enough results
                if ( results.size() == 0 )
                    {
                        System.out.printf("[%3d] rejected (all minimization failed) on microiteration %d\n", serverID, microiteration);
                        continue;
                    }
                
                // add all the results to the list of best results if applicable
                String energyString = "[";
                for (Peptide p : results)
                    energyString += String.format("%.2f, ", p.energyBreakdown.totalEnergy);
                energyString = energyString.substring(0, energyString.length()-2) + "]";
                System.out.printf("[%3d] %d structures generated on microiteration %d %s\n", serverID, results.size(), microiteration, energyString);
                
                // return the best result
                Collections.sort(results);
                return results.get(0);
        }
    }

    /** Writes the job to disk. */
    public void checkpoint()
    {
        try
            {
                 FileOutputStream fileOut = new FileOutputStream(checkpointFilename);
                 ObjectOutputStream out = new ObjectOutputStream(fileOut);
                 out.writeObject(this);
                 out.close();
                 fileOut.close();
            }
        catch (Exception e)
            {
                e.printStackTrace();
            }
    }

    @Override
    public MonteCarloResult call()
    {
        // this is the current state of the Monte Carlo simulation
        Peptide currentPeptide = startingPeptide;
        
        // perform the Metropolis Monte Carlo
        double currentAlpha = 0.0;
        Double lastBestEnergy = null;
        for (int i=0; i < maxIterations; i++)
            {
                // abort if kill file found
                if ( new File("kill.txt").isFile() )
                    {
                        System.out.println("Kill file found.  Aborting.");
                        break;
                    }

                // make a mutation (minimization is done inside this)
                Peptide candidate = mutate(currentPeptide);

                // apply modified Metropolis criterion
                boolean isAccepted = MonteCarloJob.acceptChange(currentPeptide, candidate, currentAlpha);
                double bestEnergy = bestPeptides.getBestEnergy();
                if ( i==0 )
                    lastBestEnergy = bestEnergy;
                if ( isAccepted )
                    {
                        currentPeptide = candidate;
                        if ( bestEnergy < lastBestEnergy )
                            System.out.printf("[%3d] Iter %d of %d (new best, alpha = %.6f, E = %.2f, bestE = %.2f)\n", serverID, i+1, maxIterations, currentAlpha, candidate.energyBreakdown.totalEnergy, bestEnergy);
                        else
                            System.out.printf("[%3d] Iter %d of %d (accepted, alpha = %.6f, E = %.2f, bestE = %.2f)\n", serverID, i+1, maxIterations, currentAlpha, candidate.energyBreakdown.totalEnergy, bestEnergy);
                    }
                else
                    System.out.printf("[%3d] Iter %d of %d (rejected, alpha = %.6f, E = %.2f, bestE = %.2f)\n", serverID, i+1, maxIterations, currentAlpha, candidate.energyBreakdown.totalEnergy, bestEnergy);

                // update alpha
                currentAlpha = currentAlpha + deltaAlpha;
            }
        
        // save this result to disk if requested
        if ( checkpointFilename != null )
            checkpoint();
        System.out.printf("[%3d] Finished. (bestE = %.2f)\n", serverID, bestPeptides.getBestEnergy());

        // return result
        return new MonteCarloResult(bestPeptides.getList());
    }
}
