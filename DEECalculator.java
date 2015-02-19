import java.io.*;
import java.util.*;
import com.google.common.collect.*;
import java.util.concurrent.*;
import java.util.concurrent.atomic.*;
import org.apache.commons.math3.geometry.euclidean.threed.Vector3D;
import org.apache.commons.io.FileUtils;
import org.apache.commons.math3.distribution.NormalDistribution;

/**
 * This class collects together some methods for calling DEE on a given peptide.
 */
public class DEECalculator implements Immutable
{
    /**
     * Calculates one- and two-center energies, performs DEE, and then enumerates the best poses using parallel A*.
     * Could result in no solution (causes an exception).  The best poses from A* will be returned.
     * @param peptide the peptide to pack
     * @param rotamerSpace the rotamer space to use
     * @param maxPoses the maximum number of peptides to return
     * @return the best poses, along with some randomly mutated versions of them
     */
    public static List<Peptide> generatePoses(Peptide peptide, RotamerSpace rotamerSpace, int maxPoses)
    {
        // calculate the required energies
        DEEenergyCalculator calculator = DEEenergyCalculator.analyze(rotamerSpace);
    
        // DEE
        int rounds = 0;
        int eliminatedThisRound = 0;
        int totalEliminated = 0;
        int totalSize = 0;
        DEESingles thisRound = null;
        List<List<Rotamer>> currentRotamerSpace = rotamerSpace.rotamerSpace;

        // perform first order split singles
        System.out.printf("Split singles...");
        do
            {
                if ( RotamerSpace.countRotamers(currentRotamerSpace) < 100 ) break;
                rounds++;
                thisRound = new DEESingles(calculator.rotamerSelfEnergies, calculator.rotamerInteractionEnergies,
                                           currentRotamerSpace, rotamerSpace.incompatiblePairs);
                List<List<Rotamer>> newRotamerSpace = thisRound.eliminate(DEESingles.JobType.SPLIT_FIRST_ORDER);

                // figure out how many rotamers were eliminated in this round
                eliminatedThisRound = 0;
                totalSize = 0;
                for (int i=0; i < currentRotamerSpace.size(); i++)
                    {
                        int originalSize = currentRotamerSpace.get(i).size();
                        int newSize = newRotamerSpace.get(i).size();
                        totalSize += newSize;
                        eliminatedThisRound += originalSize - newSize;
                        //System.out.printf("   %3d rotamers eliminated from position %d\n", originalSize-newSize, i);
                    }

                currentRotamerSpace = newRotamerSpace;
                totalEliminated += eliminatedThisRound;
                System.out.printf("%d, ", eliminatedThisRound);
            }
        while ( eliminatedThisRound > 0 );
        System.out.printf("done.  Eliminated %d rotamers, %d remain.\n", totalEliminated, totalSize);

        // perform magic bullet split singles
        rounds = 0;
        eliminatedThisRound = 0;
        totalEliminated = 0;
        System.out.printf("Magic bullet split singles...");
        do
            {
                if ( RotamerSpace.countRotamers(currentRotamerSpace) < 100 ) break;
                rounds++;
                thisRound = new DEESingles(calculator.rotamerSelfEnergies, calculator.rotamerInteractionEnergies,
                                                      currentRotamerSpace, rotamerSpace.incompatiblePairs);
                List<List<Rotamer>> newRotamerSpace = thisRound.eliminate(DEESingles.JobType.MB_SPLIT_FIRST_ORDER);

                // figure out how many rotamers were eliminated in this round
                eliminatedThisRound = 0;
                totalSize = 0;
                for (int i=0; i < currentRotamerSpace.size(); i++)
                    {
                        int originalSize = currentRotamerSpace.get(i).size();
                        int newSize = newRotamerSpace.get(i).size();
                        totalSize += newSize;
                        eliminatedThisRound += originalSize - newSize;
                        //System.out.printf("   %3d rotamers eliminated from position %d\n", originalSize-newSize, i);
                    }
                currentRotamerSpace = newRotamerSpace;
                totalEliminated += eliminatedThisRound;
                System.out.printf("%d, ", eliminatedThisRound);
            }
        while ( eliminatedThisRound > 0 );
        System.out.printf("done.  Eliminated %d rotamers, %d remain.\n", totalEliminated, totalSize);

        // perform magic bullet doubles
        thisRound = new DEESingles(calculator.rotamerSelfEnergies, calculator.rotamerInteractionEnergies,
                                   currentRotamerSpace, rotamerSpace.incompatiblePairs);

        rounds = 0;
        totalEliminated = 0;
        System.out.printf("Magic bullet doubles...");
        do
            {
                if ( RotamerSpace.countRotamers(currentRotamerSpace) < 100 ) break;
                rounds++;
                eliminatedThisRound = 0;
                DEEDoubles doubles = DEEDoubles.createDEEDoubles(thisRound);
                doubles.performMagicBulletDoubles();

                int oldIncompatibleSize = thisRound.incompatiblePairs.size();
                thisRound = doubles.getDEESingles();
                eliminatedThisRound = thisRound.incompatiblePairs.size() - oldIncompatibleSize;
                totalEliminated += eliminatedThisRound;
                System.out.printf("%d, ", eliminatedThisRound);
            }
        while (eliminatedThisRound > 0 && rounds <= 2 );
        System.out.printf("done.  Eliminated %d rotamer pairs, %d total.\n", totalEliminated, thisRound.incompatiblePairs.size());

        // split singles again
        rounds = 0;
        totalEliminated = 0;
        System.out.printf("Split singles repeat...");
        do
            {
                if ( RotamerSpace.countRotamers(currentRotamerSpace) < 100 ) break;
                rounds++;
                List<List<Rotamer>> newRotamerSpace = thisRound.eliminate(DEESingles.JobType.SPLIT_FIRST_ORDER);

                // figure out how many rotamers were eliminated in this round
                eliminatedThisRound = 0;
                totalSize = 0;
                for (int i=0; i < currentRotamerSpace.size(); i++)
                    {
                        int originalSize = currentRotamerSpace.get(i).size();
                        int newSize = newRotamerSpace.get(i).size();
                        totalSize += newSize;
                        eliminatedThisRound += originalSize - newSize;
                        //System.out.printf("   %3d rotamers eliminated from position %d\n", originalSize-newSize, i);
                    }
                currentRotamerSpace = newRotamerSpace;
                totalEliminated += eliminatedThisRound;
                System.out.printf("%d, ", eliminatedThisRound);
            }
        while ( eliminatedThisRound > 0 );
        System.out.printf("done.  Eliminated %d rotamers, %d remain.\n", totalEliminated, totalSize);

        // split magic double singles again
        rounds = 0;
        eliminatedThisRound = 0;
        totalEliminated = 0;
        System.out.printf("Magic bullet split singles repeat...");
        do
            {
                if ( RotamerSpace.countRotamers(currentRotamerSpace) < 100 ) break;
                rounds++;
                thisRound = new DEESingles(thisRound.rotamerSelfEnergies, thisRound.rotamerInteractionEnergies,
                                           currentRotamerSpace, thisRound.incompatiblePairs);
                List<List<Rotamer>> newRotamerSpace = thisRound.eliminate(DEESingles.JobType.MB_SPLIT_FIRST_ORDER);

                // figure out how many rotamers were eliminated in this round
                eliminatedThisRound = 0;
                totalSize = 0;
                for (int i=0; i < currentRotamerSpace.size(); i++)
                    {
                        int originalSize = currentRotamerSpace.get(i).size();
                        int newSize = newRotamerSpace.get(i).size();
                        totalSize += newSize;
                        eliminatedThisRound += originalSize - newSize;
                        //System.out.printf("   %3d rotamers eliminated from position %d\n", originalSize-newSize, i);
                    }
                currentRotamerSpace = newRotamerSpace;
                totalEliminated += eliminatedThisRound;
                System.out.printf("%d, ", eliminatedThisRound);
            }
        while ( eliminatedThisRound > 0 );
        System.out.printf("done.  Eliminated %d rotamers, %d remain.\n", totalEliminated, totalSize);

        // check there is a possible solution
        RotamerSpace.checkRotamerSpace(currentRotamerSpace, peptide, thisRound.incompatiblePairs);

        // print out the rotamer space after elimination
        System.out.println("Rotamer space after elimination: ");
        for (int i=0; i < currentRotamerSpace.size(); i++)
            {
                List<Rotamer> list = currentRotamerSpace.get(i);
                System.out.printf("   [ position %2d, %2d rotamers ]   ", i, list.size());
                for (Rotamer r : list)
                    {
                        String description = r.description;
                        String[] fields = description.split("_");
                        System.out.printf("%s ", fields[fields.length-1]);
                    }
                System.out.println();
            }
        /*totalSize = 0;
        int maxSize = currentRotamerSpace.get(0).size();
        rotamerSpace:
        for (int i=0; i < currentRotamerSpace.size(); i++)
            {
                List<Rotamer> list = currentRotamerSpace.get(i);
                totalSize += list.size();
                if ( list.size() > maxSize )
                    maxSize = list.size();
                System.out.printf("   [ position %2d, %2d rotamers ]   ", i, list.size());
                Map<AminoAcid,Boolean> map = new LinkedHashMap<>();
                aa:
                for (AminoAcid aa : AminoAcid.values())
                    {
                        if ( aa == AminoAcid.DPRO ||
                             aa == AminoAcid.LYS || aa == AminoAcid.MET || aa == AminoAcid.CYS )
                            continue;
                        boolean found = false;
                        for (Rotamer r : list)
                            {
                                if ( r.description.indexOf("transition_state") > -1 )
                                    continue;
                                else if ( r.protoAminoAcid.r.aminoAcid == aa )
                                    found = true;
                            }
                        map.put(aa,found);
                    }
                for ( Rotamer r : list )
                    {
                        if ( r.protoAminoAcid.r.description.indexOf("transition_state") > -1 )
                            {
                                System.out.println(r.protoAminoAcid.r.description.trim());
                                continue rotamerSpace;
                            }
                    }
                for (AminoAcid aa : map.keySet())
                    {
                        if ( map.get(aa) )
                            System.out.printf("%5s ", aa.shortName);
                        else
                            System.out.printf("%5s ", "");
                    }
                System.out.println();
            }*/
        RotamerSpace.printRotamerSizes(currentRotamerSpace);
        System.out.println(thisRound.incompatiblePairs.size() + " incompatible pairs.");

        // A* graph traversal
        System.out.println("\nA* Search");
        RotamerIterator iterator = new RotamerIterator(currentRotamerSpace, thisRound.rotamerSelfEnergies,
                                                       thisRound.rotamerInteractionEnergies, maxPoses);
        List<RotamerIterator.Node> solutions = iterator.iterate();
        if ( solutions.size() == 0 )
            throw new IllegalArgumentException("no valid solutions found"); 
        /*System.out.println(solutions.get(0));
        Peptide bestPeptide = Rotamer.reconstitute(peptide, solutions.get(0).rotamers);
        GaussianInputFile gjf = new GaussianInputFile(bestPeptide);
        gjf.write("test_peptides/best.gjf");
        TinkerXYZInputFile xyz = new TinkerXYZInputFile(bestPeptide);
        xyz.write("test_peptides/best.xyz");*/

        // take the best poses
        List<Future<Result>> futures = new LinkedList<>();
        List<Peptide> poses = Collections.synchronizedList(new ArrayList<Peptide>(maxPoses));
        for (RotamerIterator.Node node : solutions)
            {
                if ( futures.size() >= maxPoses )
                    break;
                List<Rotamer> rotamers = node.rotamers;
                PeptideJob job = new PeptideJob(rotamers, peptide, poses);
                Future<Result> f = GeneralThreadService.submit(job);
                futures.add(f);
            }

        // make the peptides in parallel
        GeneralThreadService.silentWaitForFutures(futures);
        System.out.printf("%d poses generated.\n", poses.size());
        return poses;
    }

    /** Reconstitutes a set of rotamers in parallel. */
    public static class PeptideJob implements WorkUnit
    {
        public final List<Rotamer> rotamers;
        public final Peptide templatePeptide;
        public final List<Peptide> targetList;

        public PeptideJob(List<Rotamer> rotamers, Peptide templatePeptide, List<Peptide> targetList)
        {
            this.rotamers = rotamers;
            this.templatePeptide = templatePeptide;
            this.targetList = targetList;
        }

        public Result call()
        {
            Peptide newPeptide = Rotamer.reconstitute(templatePeptide, rotamers);
            targetList.add(newPeptide);
            return null;
        }
    }

    /** for testing */
    public static void main(String[] args)
    {
        DatabaseLoader.go();
        List<Peptide> sheets = BetaSheetGenerator.generateSheets(5, 10, 10000, 0.01);
        sheets = new ArrayList<>(sheets);
        Collections.sort(sheets);

        // remove duplicates
        int maxResults = 64;
        double RMSDthreshold = 1.5;
        List<Peptide> results = new ArrayList<>();
        results.add(sheets.get(0));
        for (int i=1; i < sheets.size(); i++)
            {
                Peptide candidatePeptide = sheets.get(i);
                Superposition superposition = Superposition.superimpose(candidatePeptide, results);
                boolean accept = true;
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
        Peptide.writePeptideGJFs(results, "test_peptides/sheet_", 3);

        // find interesting tuples
        List<Peptide> interestingPeptides = Collections.synchronizedList(new ArrayList<Peptide>());
        List<Future<Result>> futures = new ArrayList<>();
        for (Peptide peptide : results)
            {
                RotamerFactory.InterestingJob job = new RotamerFactory.InterestingJob(peptide, interestingPeptides);
                Future<Result> f = GeneralThreadService.submit(job);
                futures.add(f);
            }
        GeneralThreadService.waitForFutures(futures);
        if ( interestingPeptides.size() == 0 )
            throw new IllegalArgumentException("no interesting peptides found");

        // write out results
        System.out.printf("\n%d peptides generated\n", interestingPeptides.size());
        Peptide.writePeptideGJFs(interestingPeptides, "test_peptides/result_", 3);
        Peptide.writePeptideCHKs(interestingPeptides, "test_peptides/result_", 3);

        Peptide peptide = interestingPeptides.get(0);
        new GaussianInputFile(peptide).write("test_peptides/interesting_peptide.gjf");
        peptide.checkpoint("test_peptides/interesting_peptide.chk");
        
        //Peptide peptide = Peptide.load("test_peptides/interesting_peptide.chk");
        peptide = HydrogenBondMutator.unmutate(peptide);
        CatalystRotamerSpace catalystRotamerSpace = new CatalystRotamerSpace(peptide,true);
        List<Peptide> poses = generatePoses(peptide, catalystRotamerSpace, 10);
        Peptide.writePeptideGJFs(poses, "test_peptides/poses_", 3);
    }
}
