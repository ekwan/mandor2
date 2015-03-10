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

    public static class LoadUnit implements WorkUnit
    {
        public final Map<BackboneFingerprint,Peptide> targetMap;
        public final String filename;

        public LoadUnit(Map<BackboneFingerprint,Peptide> targetMap, String filename)
        {
            this.targetMap = targetMap;
            this.filename = filename;
        }

        public Result call()
        {
            System.out.printf("Loading checkpoint from %s.\n", filename);
            PeptideArchive archive = PeptideArchive.load(filename);
            System.out.printf("Done loading checkpoint from %s.\n", filename);
            int numberOfPeptides = 0;
            int progress = 0;
            for (String s : archive.map.keySet())
                {
                    List<Peptide> list = archive.map.get(s);
                    numberOfPeptides += list.size();
                    ProtoAminoAcid glycine = ProtoAminoAcidDatabase.getTemplate("standard_glycine");
                    for (Peptide p : list)
                        {
                            progress++;
                            // mutate non-hairpin positions to glycine
                            Peptide p2 = p;
                            for (int i=0; i < p2.sequence.size(); i++)
                                {
                                    Residue residue = p2.sequence.get(i);
                                    if ( residue.isHairpin )
                                        continue;
                                    p2 = SidechainMutator.mutateSidechain(p2, residue, glycine);
                                }

                            BackboneFingerprint fingerprint = new BackboneFingerprint(p2);
                            targetMap.put(fingerprint,p2);
                            if ( progress % 100 == 0 )
                                System.out.printf("Processed %d of %d peptides from %s.\n", progress, numberOfPeptides, filename);
                        }
                }
            System.out.printf("Done adding unique peptides from %s.\n", filename);
            return null;
        }
    }

    /** for testing */
    public static void main(String[] args)
    {
        /*PeptideArchive archive1 = PeptideArchive.load("checkpoints/interesting.chk");
        for (String s : archive1.map.keySet())
            {
                List<Peptide> list = archive1.map.get(s);
                Peptide.writeGJFs(list, "test_peptides/designs_", 3, 1000);
                System.exit(1);
            }
        
        // make some beta sheet templates
        Map<BackboneFingerprint,Peptide> sheetMap = new ConcurrentHashMap<>(10000000, 0.75F, Settings.NUMBER_OF_THREADS);
        List<Future<Result>> futures = new ArrayList<>();
        for (File f : new File("checkpoints/").listFiles())
            {
                if ( ! f.getName().endsWith(".chk") )
                    continue;
                LoadUnit job = new LoadUnit(sheetMap, "checkpoints/" + f.getName());
                Future<Result> future = GeneralThreadService.submit(job);
                futures.add(future);
            }
        GeneralThreadService.silentWaitForFutures(futures);
        System.out.printf("Done loading.  %d unique peptides found.\n", sheetMap.size());
        */
        /*
        int numberOfSheetIterations = 20;
        int offset = 3;
        for (int i=0; i < numberOfSheetIterations; i++)
            {
                System.out.printf("=== Beta Sheet MC Run %d of %d (%d results so far) ===\n", i+1, numberOfSheetIterations, sheetMap.size());
                List<Peptide> sheets = BetaSheetGenerator.generateSheets(6, 500, 100000, 0.005);
                //List<Peptide> sheets = BetaSheetGenerator.generateSheets(5, 20, 100000, 0.0001);
                for (Peptide p : sheets)
                    {
                        BackboneFingerprint fingerprint = new BackboneFingerprint(p);
                        sheetMap.put(fingerprint,p);
                    }
                System.out.print("Serializing...");
                Map<String,List<Peptide>> backup = new HashMap<>();
                backup.put("MC run " + i, sheets);
                PeptideArchive archive = new PeptideArchive(backup);
                archive.checkpoint(String.format("checkpoints/sheet_run_%d.chk", i+offset));
                System.out.println("done.");
            }*/
        /*
        System.out.print("Creating peptide list...");
        List<Peptide> sheets = new ArrayList<>(sheetMap.values());
        System.out.println("done.");
        sheetMap = null;

        System.out.print("Garbage collecting...");
        System.gc();
        System.out.println("done.");
        */


        /*// remove duplicates
        Collections.sort(sheets);
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
        /*DatabaseLoader.go();
        List<Peptide> interestingPeptides = Collections.synchronizedList(new ArrayList<Peptide>());
        futures = new ArrayList<>();
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
        System.out.printf("\n%d interesting peptides generated\n", interestingPeptides.size());
        System.out.print("Saving interesting peptides...");
        Map<String,List<Peptide>> tempMap = new HashMap<>();
        //List<Peptide> tempList = new ArrayList<>(1000);
        //int archiveCount = 0;
        //for (Peptide p : interestingPeptides)
        //    {
        //    }
        tempMap.put("a", interestingPeptides);
        PeptideArchive peptideArchive = new PeptideArchive(tempMap);
        peptideArchive.checkpoint("checkpoints/interesting.chk");
        System.out.println("saved.");
        System.exit(1);*/

        
        //Peptide.writeGJFs(interestingPeptides, "test_peptides/result_", 3, 10);
        //Peptide.writeCHKs(interestingPeptides, "test_peptides/result_", 3, 10);
  /*      
        // mutate all non-active site positions to alanine
        DatabaseLoader.go();
        System.out.print("Loading...");
        PeptideArchive archive = PeptideArchive.load("checkpoints/interesting.chk");
        List<Peptide> interestingPeptides = new ArrayList<>();;
        for (String s : archive.map.keySet())
            {
                for (Peptide p : archive.map.get(s))
                    interestingPeptides.add(p);
            }
        System.out.println("done.");

        ProtoAminoAcid paaTemplate = ProtoAminoAcidDatabase.getTemplate("standard_alanine");
        List<Peptide> interestingPeptides2 = new ArrayList<>(interestingPeptides.size());
        for (Peptide p : interestingPeptides)
            {
                System.out.printf("%d of %d peptides mutated\r", interestingPeptides2.size(), interestingPeptides.size());
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
        System.out.println();

        // minimize poses with AMOEBA in the gas phase
        System.out.println("Minimizing all poses:");
        List<Peptide> minimizedPoses = TinkerJob.minimize(interestingPeptides2, Forcefield.AMOEBA, 3000, false, false, false, false, false);
        System.out.println("\nDone minimizing.");
        //Peptide.writeGJFs(minimizedPoses, "test_peptides/minimized_poses_", 5, maxPoses);

        // check the poses are actually sheets
        */
        /*
        PeptideArchive archive1 = PeptideArchive.load("checkpoints/designs.chk");
        List<Peptide> minimizedPoses = new ArrayList<>();
        for (String s : archive1.map.keySet())
            minimizedPoses.addAll(archive1.map.get(s));
        Map<String,List<Peptide>> resultMap = new HashMap<>(); // map from indices of TS,his,arg to peptides
        for (Peptide p : minimizedPoses)
            {
                // check to see if the structure still looks catalytic
                if ( ! isCatalytic(p) )
                   continue;
                
                // adjust peptide energy by reference energy
                //double referenceEnergy = Interaction.getAMOEBAReferenceEnergy(p);
                //EnergyBreakdown energyBreakdown = p.energyBreakdown.addReferenceEnergy(referenceEnergy);
                //Peptide adjustedPeptide = p.setEnergyBreakdown(energyBreakdown);

                // place in result map
                String signature = DEECalculator.getSignature(p);
                List<Peptide> list = resultMap.get(signature);
                if ( list == null )
                    {
                        list = new ArrayList<>();
                        resultMap.put(signature,list);
                    }
                list.add(p);
            }
        Map<Double,String> reportMap = new TreeMap<>();
        for (String s : resultMap.keySet())
            {
                List<Peptide> list = resultMap.get(s);
                Collections.sort(list);
                Peptide p = list.get(0);
                String reportString = String.format("%20s %5d  %7.2f", s, list.size(), p.energyBreakdown.totalEnergy);
                reportMap.put(p.energyBreakdown.totalEnergy, reportString);
                Peptide.writeGJFs(list, "test_peptides/design_" + s + "_", 2, 1);
                Peptide.writeCHKs(list, "test_peptides/design_" + s + "_", 2, 1);
            }
        for (Double d : reportMap.keySet())
            System.out.println(reportMap.get(d));

        //PeptideArchive archive2 = new PeptideArchive(resultMap);
        //archive2.checkpoint("checkpoints/designs.chk");
*/

        // do MC packing
        List<Future<Result>> futures = new ArrayList<>();
        for (File f : new File("test_peptides/").listFiles())
            {
                if ( !f.getName().endsWith(".chk") || f.getName().indexOf("design") == -1 )
                    continue;
                Peptide p = Peptide.load("test_peptides/" + f.getName());
                String filename = String.format("checkpoints/vsmcjob_%s", f.getName());
                System.out.println(filename);
                VariableSequenceMonteCarloJob job = null;
                if ( new File(filename).exists() )
                    {
                        job = VariableSequenceMonteCarloJob.load(filename);
                        if ( job != null )
                            {
                                if ( job.isDone() )
                                    {
                                        System.out.printf("Job from %s is complete, so doing nothing.\n", filename, job.currentIteration);
                                        continue;
                                    }
                                else
                                    System.out.printf("Resuming job from %s (%d iterations complete).\n", filename, job.currentIteration);
                            }
                    }
                if ( job == null )
                    {
                        job = new VariableSequenceMonteCarloJob(p, 0.01, 100, 100, 4, filename);
                        System.out.printf("Created new job for %s.\n", filename);
                    }
                Future<Result> f = GeneralThreadService.submit(job);
                futures.add(f);
            }
        GeneralThreadService.silentWaitForFutures(futures);
        System.out.println("All MC jobs complete.");

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
