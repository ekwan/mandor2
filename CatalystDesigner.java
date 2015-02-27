import java.io.*;
import java.util.*;
import com.google.common.collect.*;
import java.util.concurrent.*;
import java.util.concurrent.atomic.*;
import org.apache.commons.math3.geometry.euclidean.threed.Vector3D;
import org.apache.commons.io.FileUtils;
import org.apache.commons.math3.distribution.NormalDistribution;

/**
 * Tries to design beta-sheets for catalysis.  The process is:
 * 1. Make beta sheets.
 * 2. Find interesting tuples.
 * 3. Mutate non-catalytic residues to random other residues.
 * 4. Use the structure from 3 as the seed for a Monte Carlo search where
 *    mutations are single point mutation at non-catalytic residues.  Structures
 *    where the catalytic site are not preserved are rejected.  The energy is
 *    the AMOEBA energy in the gas phase + the approximate solvation energy -
 *    the reference energies for the noncatalytic residues, which should include
 *    solvation.  On each mutation, rotamer packing with A* is performed.
 */
public class CatalystDesigner implements Immutable
{
    public static boolean isCatalytic(Peptide p)
    {
        boolean isSheet = BetaSheetGenerator.isSheet(p,1);
        boolean hasBackboneContact = DEECalculator.hasBackboneContact(p);
        boolean hasArgContact = DEECalculator.hasArgContact(p);
        String signature = DEECalculator.getSignature(p);
        boolean pass = isSheet && hasBackboneContact && hasArgContact;
        System.out.printf("%10s sheet: %5b   backbone: %5b   arg: %5b  pass: %5b\n", signature, isSheet, hasBackboneContact, hasArgContact, pass);
        return pass;    
    }

    /** for testing */
    public static void main(String[] args)
    {
        // make some beta sheet templates
        DatabaseLoader.go();
        List<Peptide> sheets = BetaSheetGenerator.generateSheets(5, 50, 100000, 0.001);
        Collections.sort(sheets);

        /*// remove duplicates
        int maxResults = 10000;
        double RMSDthreshold = 0.50;
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
        */

        // find interesting tuples
        List<Peptide> interestingPeptides = Collections.synchronizedList(new ArrayList<Peptide>());
        List<Future<Result>> futures = new ArrayList<>();
        for (Peptide peptide : sheets)
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

        // mutate all non-active site positions to alanine
        ProtoAminoAcid paaTemplate = ProtoAminoAcidDatabase.getTemplate("standard_alanine");
        List<Peptide> interestingPeptides2 = new ArrayList<>(interestingPeptides.size());
        for (Peptide p : interestingPeptides)
            {
                Peptide p2 = p;
                for (int i=0; i < p.sequence.size(); i++)
                    {
                        Residue r = p2.sequence.get(i);
                        if ( r.isHairpin )
                            continue;
                        if ( r.aminoAcid == AminoAcid.GLY )
                            p2 = SidechainMutator.mutateSidechain(p2,r,paaTemplate);
                    }
                interestingPeptides2.add(p2);
            }

        // minimize poses with AMOEBA in the gas phase, then analyze with approximate solvation
        System.out.println("Minimizing all poses:");
        List<Peptide> minimizedPoses = TinkerJob.minimize(interestingPeptides2, Forcefield.AMOEBA, 2000, false, false, true, false, true);
        System.out.println("\nDone minimizing.");
        //Peptide.writeGJFs(minimizedPoses, "test_peptides/minimized_poses_", 5, maxPoses);

        // check the poses are actually sheets
        Map<String,List<Peptide>> resultMap = new HashMap<>(); // map from indices of TS,his,arg to peptides
        Peptide singlePeptide = null;
        for (Peptide p : minimizedPoses)
            {
                // check to see if the structure still looks catalytic
                if ( ! isCatalytic(p) )
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
                singlePeptide = adjustedPeptide;
                break;
            }

        // sort by design type
        
        // do MC packing
        if ( singlePeptide == null )
            throw new IllegalArgumentException("nothing to do");
        VariableSequenceMonteCarloJob job = new VariableSequenceMonteCarloJob(singlePeptide, 0.001, 10, 100, 4, "test_peptides/test.chk");
        job.call();
        Peptide.writeGJFs(job.bestPeptides.getList(), "test_peptides/vsmcjob_", 2, 100);

        // remove duplicates

/*
        // write out the results
        for (String signature : resultMap.keySet())
            {
                List<Peptide> list = resultMap.get(signature);
                Collections.sort(list); // sort by energy
                Peptide.writeGJFs(list, "test_peptides/good_" + signature + "_", 3, 1000);
                Peptide.writeCHKs(list, "test_peptides/good_" + signature + "_", 3, 1000);
            }
*/
    }
}
