import java.io.*;
import java.util.*;
import com.google.common.collect.*;
import java.util.concurrent.*;
import java.util.concurrent.atomic.*;
import org.apache.commons.math3.geometry.euclidean.threed.Vector3D;
import org.apache.commons.io.FileUtils;
import org.apache.commons.math3.distribution.NormalDistribution;

public class FixedSequenceMonteCarloJob extends MonteCarloJob
{
    /** The number of lowest energy peptides to keep */ 
    public final int maxSize;

    /** If a non-hairpin proline is selected for mutation, the chance of rotating it to the cis conformation. */
    public static final double CIS_PROLINE_PROBABILITY = 0.10;

    /** For setting omegas. */
    public static final NormalDistribution TRANS_DISTRIBUTION = new NormalDistribution(180.0,5.0);
    
    /** For setting omegas. */
    public static final NormalDistribution CIS_DISTRIBUTION = new NormalDistribution(0.0,5.0);

    /** A list of backbone atom pairs to check for clashes. */
    public static final List<Pair<Integer, Integer>> backbonePairs;
    
    /**
     * The allowed (phi,psi) angles for the hairpin D-proline and glycine.  Entries 0 and 1 are phi and psi for
     * the proline, respectively.  Entries 2 and 3 are phi and psi for the glycine, respectively.
     */
    public static final DiscreteProbabilityDistribution<List<Double>> HAIRPIN_DISTRIBUTION;

    public FixedSequenceMonteCarloJob(Peptide startingPeptide, double deltaAlpha, int maxIterations, int maxSize)
    {
        super(startingPeptide, deltaAlpha, maxIterations);
        this.maxSize = maxSize;
        this.backbonePairs = startingPeptide.getBackbonePairs();

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
        
        // choose a residue at random
        int randomIndex = -1;
        Residue residue = null;

        while (true)
            {
                randomIndex = random.nextInt(peptide.sequence.size());
                residue = peptide.sequence.get(randomIndex);
                
                // if this is a hairpin residue, give a high chance of rolling again
                if ( residue.description.indexOf("hairpin") > -1 && random.nextDouble() > 0.05 )
                    continue;
                break;
            }

        // perturb the backbone
        if ( residue.description.indexOf("hairpin") > -1 )
            {
                // if this is a hairpin residue, set the hairpin
                // this is hard-coded for a 12-mer
                
                // choose the phi,psi angles
                List<Double> hairpinAngles = HAIRPIN_DISTRIBUTION.getRandom();
                double prolinePhiAngle = hairpinAngles.get(0);
                double prolinePsiAngle = hairpinAngles.get(1);
                double glycinePhiAngle = hairpinAngles.get(2);
                double glycinePsiAngle = hairpinAngles.get(3);

                System.out.printf("[ %d ] Set hairpin backbone to: %.0f, %.0f, %.0f, %.0f.\n", serverID, prolinePhiAngle, prolinePsiAngle, glycinePhiAngle, glycinePsiAngle);

                // choose the omega angles
                double prolineOmegaAngle = BetaSheetGenerator.angleModulus(TRANS_DISTRIBUTION.sample());
                double glycineOmegaAngle = BetaSheetGenerator.angleModulus(TRANS_DISTRIBUTION.sample());

                // apply the changes
                Residue proline = newPeptide.sequence.get(5);
                newPeptide = BackboneMutator.mutatePhiPsi(newPeptide, proline, prolinePhiAngle, prolinePsiAngle);
                proline = newPeptide.sequence.get(5);
                newPeptide = BackboneMutator.mutateOmega(newPeptide, proline, prolineOmegaAngle);

                Residue glycine = newPeptide.sequence.get(6);
                newPeptide = BackboneMutator.mutatePhiPsi(newPeptide, glycine, glycinePhiAngle, glycinePsiAngle);
                glycine = newPeptide.sequence.get(6);
                newPeptide = BackboneMutator.mutateOmega(newPeptide, glycine, glycineOmegaAngle);
                
                // return the result
                return newPeptide;
            }
        else
            {
                // otherwise, choose a random (omega,phi,psi)
                System.out.printf("[ %d ] Setting backbone angles for residue %d (%s).\n", serverID, randomIndex, residue.aminoAcid.shortName);
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
                        newPeptide = BackboneMutator.mutateOmega(newPeptide, residue, omegaAngle);

                        // choose a random (phi,psi)
                        newPeptide = BackboneMutator.mutatePhiPsi(newPeptide,randomIndex+1);
                    
                        // return the result
                        return newPeptide;
                    }
                else
                    {
                        // choose a random (phi,psi)
                        newPeptide = BackboneMutator.mutateOmega(newPeptide,residue);
                        residue = newPeptide.sequence.get(randomIndex);
                        newPeptide = BackboneMutator.mutatePhiPsi(newPeptide,residue);
                        return newPeptide;
                    }
        }
            
        // Rotamer Packing

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
                // Return first peptide in solution set which is the lowest energy peptide
                return thisPeptide;
            }

    }

    /** 
    * Performs OPLS minimization with single point solvation energy.
    * @param peptide the input peptide that will be minimized
    * @return OPLS-minimized peptide
    */
    public Peptide minimizeSingle(Peptide peptide)
    {
        // Minimize each peptide on OPLS forcefield
        TinkerJob job = new TinkerJob(peptide, Forcefield.OPLS, 2000, false, true, false, false, false);
        TinkerResult result = job.call();
        return result.minimizedPeptide;
    }

    public MonteCarloResult call()
    {
        List<Peptide> bestResults = minimize(maxSize);
        return MonteCarloResult(bestResults);
    }
       
}
