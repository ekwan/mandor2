import java.util.*;
import java.io.*;
import java.util.concurrent.atomic.*;
import org.apache.commons.math3.geometry.euclidean.threed.*;

/**
 * This class represents a TINKER analysis job.
 * It is assumed that TINKER is in the path.
 * Note that this will not work with "solvate gk" because of a bug in Tinker.
 * It will work for "solvate gb," but because the tinker solvation terms are arbitrarily
 * partitioned, using this might have weird consequences.
 */
public class TinkerAnalysisJob implements WorkUnit, Serializable, Immutable
{
    /** for serialization */
    public static final long serialVersionUID = 1L;

    /** counter for generating unique filename */
    public static final AtomicInteger index = new AtomicInteger();

    /** default keywords for AMOEBA */
    public static final String DEFAULT_AMOEBA_KEYWORDS = "parameters amoebapro13.prm\n\n";

    /** default keywords for OPLS */
    public static final String DEFAULT_OPLS_KEYWORDS = "parameters oplsaal.prm\n\n";

    /** the xyz file that tinker will read from */
    public final TinkerXYZInputFile tinkerXYZInputFile;

    /** the key file that tinker will read from */
    public final TinkerKeyFile tinkerKeyFile;

    /** the peptide for which the Analysis job is being carried out */
    public final Peptide peptide;

    /**
     * Creates a job for minimizing some peptide.
     * @param peptide the peptide whose energy is to be analyzed
     * @param forcefield which forcefield to minimize on
     * @param extraKeywords any extra keywords that are desired (newlines required)
     */
    public TinkerAnalysisJob(Peptide peptide, Forcefield forcefield, String extraKeywords)
    {
        if ( peptide == null )
            throw new NullPointerException("molecule cannot be null when constructing a tinker minimization job");
        this.peptide = peptide;
        String keywords = "";
        if ( forcefield == Forcefield.AMOEBA )
            {
                keywords = DEFAULT_AMOEBA_KEYWORDS + extraKeywords;
       	        tinkerXYZInputFile = new TinkerXYZInputFile(peptide, Forcefield.AMOEBA);
            }
        else if ( forcefield == Forcefield.OPLS )
            {
                keywords = DEFAULT_OPLS_KEYWORDS + extraKeywords;
                tinkerXYZInputFile = new TinkerXYZInputFile(peptide, Forcefield.OPLS);
            }
       else
            throw new IllegalArgumentException("unrecognized forcefield type");
	    tinkerKeyFile = new TinkerKeyFile(keywords);
    }
    
    /** will run with standard keywords only */
    public TinkerAnalysisJob(Peptide peptide, Forcefield forcefield)
    {
        this(peptide,forcefield,"");
    }

    /**
     * Auto-selects a filename and runs the minimization calculation.
     * It is assumed that the TINKER binaries are in the path.
     * @return the result of the calculation
     */
    public TinkerAnalysisResult call()
    {
        // choose an appropriate base filename
        // try up to TINKER_ANALYSIS_MAX_FILENAMES times to make a unique set of filenames
        String baseFilename = "";
        counting:
        for (int i=0; i < Settings.TINKER_ANALYSIS_MAX_FILENAMES; i++)
            {
                // get a new ID number for this job
                int currentIndex = index.getAndIncrement();
                baseFilename = String.format("%s_tinker_analysis_job_%010d", Settings.HOSTNAME, currentIndex);

                // don't allow this choice of filenames if any files with this prefix already exist
                for ( File f : new File(Settings.TINKER_ANALYSIS_JOB_DIRECTORY).listFiles() )
                    {
                        if ( f.getName().startsWith(baseFilename) )
                            {
                                baseFilename = "";
                                continue counting;
                            }
                    }
                break;
            }
        if ( baseFilename.length() == 0 )
            throw new IllegalArgumentException("Unable to set filename!");

        // reset counter if necessary
        if ( index.get() > Settings.TINKER_ANALYSIS_MAX_FILENAMES )
            index.getAndSet(0);

	    // write input files to disk
        tinkerXYZInputFile.write( Settings.TINKER_ANALYSIS_JOB_DIRECTORY + baseFilename + ".xyz" );
	    tinkerKeyFile.write     ( Settings.TINKER_ANALYSIS_JOB_DIRECTORY + baseFilename + ".key" );

        // call minimize
        double elapsedTime = 0.0;
        int exitValue = -1;
        //boolean badGeometry = false;
        try
            {
                if ( Settings.PLATFORM == Settings.Platform.DOS )
                    throw new IllegalArgumentException("mandor doesn't work on DOS yet");
                else if ( Settings.PLATFORM == Settings.Platform.LINUX )
                    {
                        long startTime = System.currentTimeMillis();
                        String runString = Settings.TINKER_ANALYSIS_JOB_DIRECTORY + "run_tinker_analysis.sh " +
                                           Settings.TINKER_ANALYSIS_JOB_DIRECTORY + " " + baseFilename;
                        //System.out.println(runString);
                        Process process = Runtime.getRuntime().exec(runString);
                        process.waitFor();
                        long endTime = System.currentTimeMillis();
                        exitValue = process.exitValue();
                        elapsedTime = (endTime - startTime) / 1000.0;
                    }
            }
        catch (Exception e)
            {
                System.out.println("Error while running Tinker job:");
                System.out.println(baseFilename);
                e.printStackTrace();
            }

        // remove files on abnormal termination
        if ( exitValue != 0 )
            {
                try
                    {
                        File[] files = new File(Settings.TINKER_ANALYSIS_JOB_DIRECTORY).listFiles();
                        for ( File f : files )
                            {
                                String filename = f.getName();
                                if ( f.getName().startsWith(baseFilename) )
                                    {
                                        //f.delete();
                                        //System.out.println(f.getName() + " deleted.");
                                    }
                            }
                    }
                catch (Exception e)
                    {
                        System.out.println("Error while trying to delete files:");
                        e.printStackTrace();
                    }
            }

        // check if the job completed correctly
        if ( exitValue == -1 )
            throw new IllegalArgumentException(baseFilename + " exceeded the allotted time");
        //else if ( badGeometry )
        //    throw new IllegalArgumentException(baseFilename + ": interpolation error");
        else if ( exitValue != 0 )
            {
                String tail = "";
                try
                    {
                        OutputFileFormat errorOutput = new OutputFileFormat(Settings.TINKER_ANALYSIS_JOB_DIRECTORY + baseFilename + ".out") {};
                        String[] lines = errorOutput.stringRepresentation.split("\n");
                        int length = lines.length;
                        for (int i=Math.max(length-10,0); i < length; i++)
                            tail += lines[i] + "\n";
                    }
                catch (Exception e)
                    {
                    }
                throw new IllegalArgumentException("error code " + exitValue + " while runing tinker analysis job: " + baseFilename + "!\n" + tail);
            }

        // retrieve output analysis file
        TinkerAnalyzeOutputFile analysisOutput = new TinkerAnalyzeOutputFile(Settings.TINKER_ANALYSIS_JOB_DIRECTORY + baseFilename + ".txt", peptide);
        
        // delete files after normal termination
        try
            {
                File[] files = new File(Settings.TINKER_ANALYSIS_JOB_DIRECTORY).listFiles();
                for ( File f : files )
                    {
                        String filename = f.getName();
                        if ( f.getName().startsWith(baseFilename) )
                            {
                                f.delete();
                                //System.out.println(f.getName() + " deleted.");
                            }
                    }
            }
        catch (Exception e)
            {
                System.out.println("Error while trying to delete files:");
                e.printStackTrace();
            }

        // construct and return result
        return new TinkerAnalysisResult(analysisOutput, elapsedTime);
    }
    
    /**
     * A class representing output files from a analysis call.
     * Contains the TinkerAnalyzeOutputFile produced from the analysis call.
     */
    public static class TinkerAnalysisResult implements Result, Serializable
    {
        /** for serialization */
        public static final long serialVersionUID = 1L;

        /** The Tinker Analysis File that is produced by analyze containing the energy breakdown of a peptide */
        public final TinkerAnalyzeOutputFile tinkerAnalysisFile;

        /** The run time of the call in seconds. */
        public final double elapsedTime;

        public TinkerAnalysisResult(TinkerAnalyzeOutputFile tinkerAnalysisFile, double elapsedTime)
        {
            this.elapsedTime = elapsedTime;
            this.tinkerAnalysisFile = tinkerAnalysisFile;
        }
	
        @Override
        public int hashCode()
        {
            return Objects.hash(tinkerAnalysisFile, elapsedTime);
        }

        @Override
        public boolean equals(Object obj)
        {
            if ( obj == null )
                return false;
            if ( obj == this )
                return true;
            if ( !(obj instanceof TinkerAnalysisResult) )
                return false;

            TinkerAnalysisResult result = (TinkerAnalysisResult)obj;
            if ( Objects.equals(this.tinkerAnalysisFile, result.tinkerAnalysisFile) &&
                                this.elapsedTime == result.elapsedTime )
                return true;
            return false;
        }

        @Override
        public String toString()
        {
            return tinkerAnalysisFile.toString() + String.format(" time = %.3f s\n\n", elapsedTime);
        }
    }
    
    @Override
    public int hashCode()
    {
        return Objects.hash(tinkerXYZInputFile, tinkerKeyFile, peptide);
    }

    @Override
    public boolean equals(Object obj)
    {
        if ( obj == null )
            return false;
        if ( obj == this )
            return true;
        if ( !(obj instanceof TinkerAnalysisJob) )
            return false;

        TinkerAnalysisJob job = (TinkerAnalysisJob)obj;
        if ( Objects.equals(this.tinkerXYZInputFile, job.tinkerXYZInputFile) &&
             Objects.equals(this.tinkerKeyFile,      job.tinkerKeyFile) && Objects.equals(this.peptide, job.peptide))
            return true;
        return false;
    }

    @Override
    public String toString()
    {
        return "TinkerJob:\n" + tinkerXYZInputFile.stringRepresentation + "\n" + tinkerKeyFile.stringRepresentation + "\n" + peptide.toString();
    }
  
    /** For testing. */
    public static void main(String[] args)
    {
        DatabaseLoader.go();
        System.out.println("building");
        List<ProtoAminoAcid> sequence = ProtoAminoAcidDatabase.getSpecificSequence("arg","met","standard_ala","gly","d_proline", "gly", "phe", "val", "hd", "l_pro");
        Peptide peptide = PeptideFactory.createPeptide(sequence);
        
        for (int i=0; i < peptide.sequence.size(); i++)
            {
                peptide = BackboneMutator.mutateOmega(peptide, i);
                peptide = BackboneMutator.mutatePhiPsi(peptide, i);
                peptide = RotamerMutator.mutateChis(peptide, i);
            }
        peptide = PeptideFactory.setHairpinAngles(peptide);

        System.out.println("minimizing");
        TinkerMinimizationJob job = new TinkerMinimizationJob(peptide, Forcefield.OPLS, "maxiter 2000\n\n");
        //TinkerMinimizationJob job = new TinkerMinimizationJob(peptide, Forcefield.AMOEBA, "maxiter 2000\n\n");
        TinkerMinimizationJob.TinkerMinimizationResult result = job.call();
        Molecule minimizedMolecule = result.tinkerXYZOutputFile.molecule;
        Peptide newPeptide = peptide.setPositions(minimizedMolecule);
        System.out.println("grad is " + result.tinkerMinimizationLogFile.gradient);
        System.out.println(result.tinkerMinimizationLogFile.iterations + " iterations performed");

        System.out.println("amoeba gas phase");
        TinkerAnalysisJob job2 = new TinkerAnalysisJob(newPeptide, Forcefield.AMOEBA, "\n\n");
        TinkerAnalysisJob.TinkerAnalysisResult result2 = job2.call();
        TinkerAnalyzeOutputFile outputFile2 = result2.tinkerAnalysisFile;
        System.out.println(outputFile2.energyByResidue);
        System.out.println(outputFile2.totalEnergy);
        
        System.out.println("opls with gb");
        TinkerAnalysisJob job3 = new TinkerAnalysisJob(newPeptide, Forcefield.OPLS, "solvate gb\n\n");
        TinkerAnalysisJob.TinkerAnalysisResult result3 = job3.call();
        TinkerAnalyzeOutputFile outputFile3 = result3.tinkerAnalysisFile;
        System.out.println(outputFile3.energyByResidue);
        System.out.println(outputFile3.totalEnergy);

        System.out.println("amoeba with gk");
        TinkerAnalysisJob job4 = new TinkerAnalysisJob(newPeptide, Forcefield.AMOEBA, "solvate gk\n\n");
        TinkerAnalysisJob.TinkerAnalysisResult result4 = job4.call();
        TinkerAnalyzeOutputFile outputFile4 = result4.tinkerAnalysisFile;
        System.out.println(outputFile4.energyByResidue);
        System.out.println(outputFile4.totalEnergy);
    }
}
