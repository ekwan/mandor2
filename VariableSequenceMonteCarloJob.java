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
     * @param peptide the peptide whose backbone we should perturb
     * @return the perturbed peptide
     */
    public Peptide mutate(Peptide peptide)
    {
        // random number generator
        ThreadLocalRandom random = ThreadLocalRandom.current();

        // figure out which positions can be mutated
        List<Integer> validPositions = new ArrayList<>(peptide.sequence.size());
        for (int i=0; i < peptide.sequence.size(); i++)
            {
                Residue residue = peptide.sequence.get(i);
                if ( residue.isHairpin || residue.aminoAcid == AminoAcid.TS ||
                     residue.aminoAcid == AminoAcid.HIS || residue.aminoAcid == AminoAcid.ARG )i 
                    continue;
                validPositions.add(i);
            }
 
        // mutate until valid A* poses can be generated
        int microiteration = 0;
        while (true)
            {
                microiteration++;
                FixedSequenceRotamerSpace fixedSequenceRotamerSpace = null;
                
                // choose a random position
                int randomIndex = validPositions.get(random.nextInt(validPositions.size()));
                Residue residue = peptide.sequence.get(randomIndex);

                // choose an amino acid to mutate to, it's always a different one than the one we have now
                AminoAcid currentAminoAcid = residue.aminoAcid;
                List<AminoAcid> mutationOutcomes = new ArrayList<>(CatalystRotamerSpace.MUTATION_OUTCOMES);
                mutationOutcomes.remove(currentAminoAcid);
                AminoAcid randomAminoAcid = mutationOutcomes.get(random.nextInt(mutationOutcomes.size()));
                
                // mutate at random
                int paaIndex = ProtoAminoAcidDatabase.KEYS.indexOf(randomAminoAcid);
                if ( paaIndex == -1 )
                    throw new IllegalArgumentException("no protoaminoacid templates found for " + randomAminoAcid.shortName);
                List<ProtoAminoAcid> paaList = ProtoAminoAcidDatabase.VALUES.get(paaIndex);
                ProtoAminoAcid protoAminoAcid = paaList.get(random.nextInt(paaList.size()));
                Peptide newPeptide = SidechainMutator.mutateSidechain(peptide, residue, protoAminoAcid);

                // compute the reference energy for this peptide
                double referenceEnergy = Interaction.getAMOEBAReferenceEnergy(newPeptide);

                // rotamer pack with an A* iteration
                //
                // we must call this in a try-catch clause because it's possible the peptide could be in a bad conformation
                // in which it is impossible to place any rotamers
                try
                    {
                        // note that the includeHN should be set to false if we want to do A* here
                        varaibleSequenceRotamerSpace = new VariableSequenceSequenceRotamerSpace(newPeptide, false, false);
                    }
                catch (Exception e)
                    {
                        if ( e instanceof IllegalArgumentException )
                            System.out.printf("[%3d] rejected (packing failure: %s) on microiteration %d\n", serverID, e.getMessage(), microiteration);
                        else
                            {
                                System.out.printf("[%3d] rejected (rotamer packing unknown error) on microiteration %d\n", serverID, microiteration);
                                e.printStackTrace();
                            }
                        continue;
                    }

                // perform A* iteration
                AStarEnergyCalculator calculator = AStarEnergyCalculator.analyze(variableSequenceRotamerSpace, false);
                RotamerIterator iterator = new RotamerIterator(variableSequenceRotamerSpace.rotamerSpace, calculator.rotamerSelfEnergies, calculator.rotamerInteractionEnergies, rotamersPerIteration, false);
                List<RotamerIterator.Node> solutions = iterator.iterate();
                if ( solutions.size() == 0 )
                    {
                        System.out.printf("[%3d] rejected (no A* solutions) on microiteration %d\n", serverID, microiteration);
                        continue;
                    }

                // minimize peptides in serial
                List<Peptide> results = new ArrayList<>(rotamersPerIteration+100);
                for (RotamerIterator.Node node : solutions)
                    {
                        if ( results.size() >= rotamersPerIteration )
                            break;

                        // solutions come out lowest energy first
                        List<Rotamer> rotamers = node.rotamers;
                        Peptide thisPeptide = Rotamer.reconstitute(newPeptide, rotamers);

                        // minimize with gas phase AMOEBA with approximate solvation, throw out peptides that error
                        Peptide minimizedPeptide = null;
                        try { minimizedPeptide = TinkerJob.minimize(thisPeptide, Forcefield.AMOEBA, 2000, false, false, true, true, false); }
                        catch (Exception e) { continue; }

                        // adjust total energy for reference energy
                        minimizedPeptide = minimizedPeptide.setEnergyBreakdown(energyBreakdown);         

                        results.add(minimizedPeptide);
                        bestPeptides.add(minimizedPeptide);
                   }
                //System.out.printf("[%3d] %d peptides have been minimized.\n", serverID, results.size());

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
                boolean isAccepted = MonteCarloJob.acceptChange(currentPeptide, candidate, currentAlpha);
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
}
