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
        RotamerSpace.printRotamerSizes("Rotamer space after elimination:", currentRotamerSpace);
        for (int i=0; i < currentRotamerSpace.size(); i++)
            {
                List<Rotamer> list = currentRotamerSpace.get(i);
                System.out.printf("   [ position %2d, %2d rotamers ]   ", i, list.size());
                Set<String> descriptions = new TreeSet<>();
                for (Rotamer r : list)
                    {
                        String description = r.description;
                        String[] fields = description.split("_");
                        String alteredDescription = description.replaceAll(fields[0] + "_","");
                        descriptions.add(alteredDescription);
                    }
                for (String s : descriptions)
                    System.out.print(s + " ");
                System.out.println();
            }
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
                
                // check charge
                int charge = PeptideChargeCalculator.getCharge(peptide, rotamers);
                if ( charge < 0 || charge > 1 )
                    continue;

                PeptideJob job = new PeptideJob(rotamers, peptide, poses);
                Future<Result> f = GeneralThreadService.submit(job);
                futures.add(f);
            }

        // make the peptides in parallel
        GeneralThreadService.silentWaitForFutures(futures);
        System.out.printf("%d poses generated.\n", poses.size());
        for (Peptide p : poses)
            System.out.println(p.name);
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

    /**
     * Checks to make sure that the arginine and TS oxygen are in proper contact.
     */
    public static boolean hasArgContact(Peptide peptide)
    {
        Atom TSoxygen         = RotamerFactory.locateSingleAtom(ImmutableSet.of(408,428), peptide.contents);
        List<Atom> arginineHs = RotamerFactory.locateAtoms(ImmutableSet.of(209,212), peptide.contents);
        for (Atom arginineAtom : arginineHs)
            {
                if ( Molecule.getDistance(TSoxygen, arginineAtom) < 2.2 )
                    return true;
            }
        return false;
    }

    /**
     * Ensures the transition state oxygen and a backbone HN are in proper contact.
     */
    public static boolean hasBackboneContact(Peptide peptide)
    {
        Atom TSoxygen = RotamerFactory.locateSingleAtom(ImmutableSet.of(408,428), peptide.contents);
        for (Residue r : peptide.sequence)
            {
                if ( r.HN != null && Molecule.getDistance(r.HN, TSoxygen) < 2.2 &&
                     Molecule.getAngle(r.N, r.HN, TSoxygen) > 120.0 )
                    return true;
            }
        return false;
    }

    /**
     * Returns a string of the form 0_2_4 where the numbers are the indices of the
     * TS, histidine, and arginine, respectively.  If something is not found or found twice
     * an exception is thrown.
     * @param peptide the peptide to analyze
     * @return the index string
     */
    public static String getSignature(Peptide peptide)
    {
        String returnString = "";
        int TSindex = -1;
        int histidineIndex = -1;
        int arginineIndex = -1;
        for (int i=0; i < peptide.sequence.size(); i++)
            {
                Residue residue = peptide.sequence.get(i);
                String description = residue.description;
                if ( description.indexOf("transition_state") > -1 )
                    {
                        if ( TSindex == -1 )
                            TSindex = i;
                        else throw new IllegalArgumentException("duplicate TS");
                    }
                else if ( description.indexOf("histidine") > -1 )
                    {
                        if ( histidineIndex == -1 )
                            histidineIndex = i;
                        else throw new IllegalArgumentException("duplicate histidine");
                    }
                else if ( description.indexOf("arginine") > -1 )
                    {
                        if ( arginineIndex == -1 )
                            arginineIndex = i;
                        else throw new IllegalArgumentException("duplicate arginine");
                    }
            }
        if ( TSindex == -1 || histidineIndex == -1 || arginineIndex == -1 )
            throw new IllegalArgumentException(String.format("Indices not found: TS %d, histidine %d, arginine %d", TSindex, histidineIndex, arginineIndex));
        //returnString = String.format("%d_%d_%d", TSindex, histidineIndex, arginineIndex);
        Set<Integer> set = ImmutableSet.of(TSindex, histidineIndex, arginineIndex);
        set = new TreeSet<>(set);
        for (Integer i : set)
            returnString += i + "_";
        return returnString.substring(0,returnString.length()-1);
    }

    /** for testing */
    public static void main(String[] args)
    {
        // make some beta sheet templates
        DatabaseLoader.go();
        List<Peptide> sheets = BetaSheetGenerator.generateSheets(5, 100, 100000, 0.001);
        Collections.sort(sheets);

        // remove duplicates
        int maxResults = 1000;
        double RMSDthreshold = 0.5;
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
        System.out.printf("%d unique sheets generated\n", results.size());
        //Peptide.writeGJFs(results, "test_peptides/sheet_", 3, maxResults);

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
        //Peptide.writeGJFs(interestingPeptides, "test_peptides/result_", 3, 10);
        //Peptide.writeCHKs(interestingPeptides, "test_peptides/result_", 3, 10);

        // rotamer pack
        List<Peptide> poses = new ArrayList<>();
        int maxPoses = 10000;
        for (int i=0; i < interestingPeptides.size(); i++)
            {
                if ( poses.size() > maxPoses )
                    {
                        System.out.println("Got enough poses, so not packing anymore.");
                        break;
                    }
                System.out.printf("[ Rotamer packing design %d of %d (%d poses generated so far)... ]\n", i+1, interestingPeptides.size(), poses.size());
                Peptide peptide = interestingPeptides.get(i);
                peptide = HydrogenBondMutator.unmutate(peptide);
                try
                    {
                        CatalystRotamerSpace catalystRotamerSpace = new CatalystRotamerSpace(peptide,true);
                        List<Peptide> thisPoses = generatePoses(peptide, catalystRotamerSpace, 10);
                        poses.addAll(thisPoses);
                    }
                catch (Exception e)
                    {
                        if ( e instanceof IllegalArgumentException )
                            System.out.println(e.getMessage());
                        else
                            e.printStackTrace();
                    }
            }
        //Peptide.writeGJFs(poses, "test_peptides/original_poses_", 5, maxPoses);
    
        // mutate back to close contact forcefield
        poses = HydrogenBondMutator.mutate(poses);

        // minimize poses
        System.out.println("Minimizing all poses:");
        List<Peptide> minimizedPoses = BetaSheetGenerator.minimizeSheets(poses, 2000, Forcefield.AMOEBA);
        System.out.println("\nDone minimizing.");
        //Peptide.writeGJFs(minimizedPoses, "test_peptides/minimized_poses_", 5, maxPoses);

        // check the poses are actually sheets
        Map<String,List<Peptide>> resultMap = new HashMap<>(); // map from indices of TS,his,arg to peptides
        for (Peptide p : minimizedPoses)
            {
                boolean isSheet = BetaSheetGenerator.isSheet(p,1);
                boolean hasBackboneContact = hasBackboneContact(p);
                boolean hasArgContact = hasArgContact(p);
                String signature = getSignature(p);
                boolean pass = isSheet && hasBackboneContact && hasArgContact;
                //System.out.printf("%10s sheet: %5b   backbone: %5b   arg: %5b  pass: %5b\n", signature, isSheet, hasBackboneContact, hasArgContact, pass);
                if ( !pass )
                    continue;
                
                // adjust peptide energy by reference energy
                double referenceEnergy = Interaction.getAMOEBAReferenceEnergy(p);
                EnergyBreakdown energyBreakdown = p.energyBreakdown.addReferenceEnergy(referenceEnergy);
                Peptide adjustedPeptide = p.setEnergyBreakdown(energyBreakdown);

                // place in result map
                List<Peptide> list = resultMap.get(signature);
                if ( list == null )
                    {
                        list = new ArrayList<>();
                        resultMap.put(signature,list);
                    }
                list.add(adjustedPeptide);
            }

        // write out the results
        for (String signature : resultMap.keySet())
            {
                List<Peptide> list = resultMap.get(signature);
                Collections.sort(list); // sort by energy
                Peptide.writeGJFs(list, "test_peptides/good_" + signature + "_", 3, maxPoses);
                Peptide.writeCHKs(list, "test_peptides/good_" + signature + "_", 3, maxPoses);
            }
    }
}
