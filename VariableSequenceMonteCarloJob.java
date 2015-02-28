import java.io.*;
import java.util.*;
import com.google.common.collect.*;
import java.util.concurrent.*;
import java.util.concurrent.atomic.*;
import org.apache.commons.math3.geometry.euclidean.threed.Vector3D;
import org.apache.commons.io.FileUtils;
import org.apache.commons.math3.distribution.NormalDistribution;

/**
 * Optimizes the non-active-site residues of a catalyst design.
 */
public class VariableSequenceMonteCarloJob extends MonteCarloJob implements Serializable
{
    /** For serialization. */
    public static final long serialVersionUID = 1L;

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
     * @param checkpointFilename job will be serialized here when done
     */
    public VariableSequenceMonteCarloJob(Peptide startingPeptide, double deltaAlpha, int maxIterations, int maxSize,
                                      int rotamersPerIteration, String checkpointFilename)
    {
        super(startingPeptide, deltaAlpha, maxIterations);
        this.bestPeptides = new MonteCarloJob.PeptideList(maxSize);
        this.rotamersPerIteration = rotamersPerIteration;
        this.checkpointFilename = checkpointFilename;
        this.serverID = RemoteWorkUnit.ID_GENERATOR.incrementAndGet();
    }

    /**
     * Mutates a residue at random without changing the backbone.  The procedure is:<p>
     * (1) Select a non-active-site residue at random.<p>
     * (2) Mutate to a different amino acid.<p>
     * (3) Rotamer pack, minimize, and retain only catalytic results.<p>
     * Expects and outputs close contact peptides.
     * @param peptide the peptide whose backbone we should perturb
     * @return the perturbed peptide
     */
    public Peptide mutate(Peptide inputPeptide)
    {
        Peptide.writeGJFs(bestPeptides.getList(), "test_peptides/best_", 2, 10);
        Peptide.writeCHKs(bestPeptides.getList(), "test_peptides/best_", 2, 10);

        // random number generator
        System.out.println(inputPeptide.name);
        ThreadLocalRandom random = ThreadLocalRandom.current();
        Peptide peptide = HydrogenBondMutator.unmutate(inputPeptide);

        // figure out which positions can be mutated
        List<Integer> validPositions = new ArrayList<>(peptide.sequence.size());
        for (int i=0; i < peptide.sequence.size(); i++)
            {
                Residue residue = peptide.sequence.get(i);
                if ( residue.isHairpin || residue.aminoAcid == AminoAcid.TS ||
                     residue.aminoAcid == AminoAcid.HIS || residue.aminoAcid == AminoAcid.ARG ) 
                    continue;
                validPositions.add(i);
            }
 
        // create a list of positions in the peptide sequence that can be mutated and what they can be mutated to
        List<Pair<Integer,ProtoAminoAcid>> mutationOutcomes = new ArrayList<>();
        for (Integer i : validPositions)
            {
                // make a list of all the amino acids that we can mutate to at this position
                Residue residue = peptide.sequence.get(i);
                AminoAcid currentAminoAcid = residue.aminoAcid;
                List<AminoAcid> excludeAminoAcids = new ArrayList<>();
                excludeAminoAcids.add(currentAminoAcid);
                int charge = PeptideChargeCalculator.getCharge(peptide, i);
                if ( charge == 0 )
                    {
                        excludeAminoAcids.add(AminoAcid.ASP);
                        excludeAminoAcids.add(AminoAcid.GLU);
                    }
                for (ProtoAminoAcid paa : CatalystRotamerSpace.MUTATION_OUTCOMES)
                    {
                        AminoAcid aa = paa.residue.aminoAcid;
                        if ( excludeAminoAcids.contains(aa) )
                            continue;
                        Pair<Integer,ProtoAminoAcid> pair = new Pair<>(i, paa);
                        mutationOutcomes.add(pair);
                    }
            }
        Collections.shuffle(mutationOutcomes);

        // run through the first structure in each possibility
        List<Pair<Peptide,RotamerIterator.Node>> otherSolutionsToTry = new ArrayList<>();
        for (int mutationCount=0; mutationCount < mutationOutcomes.size(); mutationCount++)
            {
                Pair<Integer,ProtoAminoAcid> mutationOutcome = mutationOutcomes.get(mutationCount);
                int i = mutationOutcome.getFirst();
                ProtoAminoAcid paa = mutationOutcome.getSecond();
                System.out.printf("Mutating position %d to %s (%d of %d possibilities).\n", i, paa.residue.description, mutationCount+1, mutationOutcomes.size());

                // make the mutation
                Residue residue = peptide.sequence.get(i);
                Peptide newPeptide = SidechainMutator.mutateSidechain(peptide, residue, paa);

                // pack the peptide
                VariableSequenceRotamerSpace variableSequenceRotamerSpace = null;
                try
                    {
                        // note that the includeHN should be set to false if we want to do A* here
                        variableSequenceRotamerSpace = new VariableSequenceRotamerSpace(newPeptide, false, false);
                    }
                catch (Exception e)
                    {
                        if ( e instanceof IllegalArgumentException )
                            System.out.printf("[%3d] microiteration rejected (packing failure: %s)\n", serverID, e.getMessage());
                        else
                            {
                                System.out.printf("[%3d] microiteration rejected (rotamer packing unknown error)\n", serverID);
                                e.printStackTrace();
                            }
                        continue;
                    }

                // check for empty nodes
                boolean empty = true;
                for (List<Rotamer> list : variableSequenceRotamerSpace.rotamerSpace)
                    {
                        if ( list.size() > 1 )
                            {
                                empty = false;
                                break;
                            }
                    }
                if ( empty )
                    {
                        System.out.printf("[%3d] microiteration rejected (empty node)\n", serverID);
                        continue;
                    }

                // rotamer pack
                AStarEnergyCalculator calculator = AStarEnergyCalculator.analyze(variableSequenceRotamerSpace, false);
                RotamerIterator iterator = new RotamerIterator(variableSequenceRotamerSpace.rotamerSpace, calculator.rotamerSelfEnergies, calculator.rotamerInteractionEnergies, rotamersPerIteration, false);
                List<RotamerIterator.Node> solutions = iterator.iterate();
                if ( solutions.size() == 0 )
                    {
                        System.out.printf("[%3d] microiteration rejected (no A* solutions)\n", serverID);
                        continue;
                    }
  
                // keep the other solutions, up to rotamersPerIteration
                List<RotamerIterator.Node> otherSolutions = new ArrayList<>();
                for (int solutionCount=1; solutionCount < Math.min(solutions.size(), rotamersPerIteration); solutionCount++)
                    otherSolutions.add(solutions.get(solutionCount));
                for (RotamerIterator.Node node : otherSolutions)
                    {
                        Pair<Peptide,RotamerIterator.Node> pair = new Pair<>(newPeptide, node);
                        otherSolutionsToTry.add(pair);
                    }
               
                // reconstitute the first pose
                RotamerIterator.Node firstNode = solutions.get(0);
                List<Rotamer> rotamers = firstNode.rotamers;
                Peptide thisPeptide = Rotamer.reconstitute(newPeptide, rotamers);
                thisPeptide = HydrogenBondMutator.mutate(thisPeptide);

                // minimize with AMOEBA with analysis (approximate OMNISOL solvation)
                Peptide minimizedPeptide = null;
                try { minimizedPeptide = TinkerJob.minimize(thisPeptide, Forcefield.AMOEBA, 250, false, false, true, false, true); }
                catch (Exception e) { e.printStackTrace(); continue; }

                // check that the geometry is still catalytic
                if ( !CatalystDesigner.isCatalytic(minimizedPeptide) )
                    continue;

                // get the reference energy
                double referenceEnergy = Interaction.getAMOEBAReferenceEnergy(minimizedPeptide);
           
                // adjust total energy for reference energy
                double referencedEnergy = minimizedPeptide.energyBreakdown.totalEnergy - referenceEnergy;
                System.out.printf("referenced energy is %.2f\n", referencedEnergy);
                EnergyBreakdown energyBreakdown = new EnergyBreakdown(null, referencedEnergy, 0.0, referencedEnergy, null, Forcefield.AMOEBA);
                minimizedPeptide = minimizedPeptide.setEnergyBreakdown(energyBreakdown);         
                bestPeptides.add(minimizedPeptide);
                return minimizedPeptide; 
            }

        // run randomly through the remaining
        Collections.shuffle(otherSolutionsToTry);
        for (int i=0; i < otherSolutionsToTry.size(); i++)
            {
                Pair<Peptide,RotamerIterator.Node> pair = otherSolutionsToTry.get(i);
                Peptide newPeptide = pair.getFirst();
                System.out.printf("Trying remaining pose %d of %d (%s).\n", i+1, otherSolutionsToTry.size(), newPeptide.name.split("@")[0]);
                RotamerIterator.Node node = pair.getSecond();
                List<Rotamer> rotamers = node.rotamers;
                Peptide thisPeptide = Rotamer.reconstitute(newPeptide, rotamers);
                thisPeptide = HydrogenBondMutator.mutate(thisPeptide);

                // minimize with AMOEBA with analysis (approximate OMNISOL solvation)
                Peptide minimizedPeptide = null;
                try { minimizedPeptide = TinkerJob.minimize(thisPeptide, Forcefield.AMOEBA, 250, false, false, true, false, true); }
                catch (Exception e) { e.printStackTrace(); continue; }

                // check that the geometry is still catalytic
                if ( !CatalystDesigner.isCatalytic(minimizedPeptide) )
                    continue;

                // get the reference energy
                double referenceEnergy = Interaction.getAMOEBAReferenceEnergy(minimizedPeptide);
           
                // adjust total energy for reference energy
                double referencedEnergy = minimizedPeptide.energyBreakdown.totalEnergy - referenceEnergy;
                System.out.printf("referenced energy is %.2f\n", referencedEnergy);
                EnergyBreakdown energyBreakdown = new EnergyBreakdown(null, referencedEnergy, 0.0, referencedEnergy, null, Forcefield.AMOEBA);
                minimizedPeptide = minimizedPeptide.setEnergyBreakdown(energyBreakdown);         
                bestPeptides.add(minimizedPeptide);
                return minimizedPeptide; 
            }

        // repack the peptide

        // we've run out of things to do, so throw an exception
        throw new IllegalArgumentException("nothing left to do");
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
        double firstEnergy = bestPeptides.getBestEnergy();
        double lastBestEnergy = firstEnergy;
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
                boolean isAccepted = true;
                if ( i > 0 ) // auto-accept first iteration
                    isAccepted = MonteCarloJob.acceptChange(currentPeptide, candidate, currentAlpha);
                double thisEnergy = candidate.energyBreakdown.totalEnergy;
                double bestEnergy = bestPeptides.getBestEnergy();
                if ( isAccepted )
                    {
                        currentPeptide = candidate;
                        if ( bestEnergy < lastBestEnergy )
                            System.out.printf("[%3d] Iter %d of %d (***new best***, alpha = %.6f, E = %.2f, bestE = %.2f, initialE = %.2f)\n", serverID, i+1, maxIterations, currentAlpha, thisEnergy, bestEnergy, firstEnergy);
                        else
                            System.out.printf("[%3d] Iter %d of %d (accepted, alpha = %.6f, E = %.2f, bestE = %.2f, initialE = %.2f)\n", serverID, i+1, maxIterations, currentAlpha, thisEnergy, bestEnergy, firstEnergy);
                    }
                else
                    System.out.printf("[%3d] Iter %d of %d (rejected, alpha = %.6f, E = %.2f, bestE = %.2f, initialE = %.2f)\n", serverID, i+1, maxIterations, currentAlpha, thisEnergy, bestEnergy, firstEnergy);

                // update alpha
                currentAlpha = currentAlpha + deltaAlpha;
                lastBestEnergy = bestEnergy;
            }
        
        // save this result to disk if requested
        if ( checkpointFilename != null )
            checkpoint();
        System.out.printf("[%3d] Finished. (bestE = %.2f)\n", serverID, bestPeptides.getBestEnergy());

        // return result
        return new MonteCarloResult(bestPeptides.getList());
    }

    public static void main(String[] args)
    {
        Peptide peptide = Peptide.load("test_peptides/singlePeptide.chk");
        VariableSequenceMonteCarloJob job = new VariableSequenceMonteCarloJob(peptide, 0.001, 1000, 100, 4, "test_peptides/test.chk");
        job.call();
        Peptide.writeGJFs(job.bestPeptides.getList(), "test_peptides/vsmcjob_", 2, 100);
    }
}
