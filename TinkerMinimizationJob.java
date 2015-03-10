import java.util.*;
import java.io.*;
import java.util.concurrent.atomic.*;
import org.apache.commons.math3.geometry.euclidean.threed.*;

/**
 * This class represents a TINKER minimization job.
 * It is assumed that TINKER is in the path.
 */
public class TinkerMinimizationJob implements WorkUnit, Serializable, Immutable
{
    /** for serialization */
    public static final long serialVersionUID = 1L;

    /** counter for generating unique filename */
    public static final AtomicInteger index = new AtomicInteger();

    /** default keywords for AMOEBA */
    public static final String DEFAULT_AMOEBA_KEYWORDS = "parameters amoebapro13.prm\nwriteout 200\n\n";

    /** default keywords for OPLS */
    public static final String DEFAULT_OPLS_KEYWORDS = "parameters oplsaal.prm\nwriteout 200\n\n";

    /** jobs are killed after this length of time in seconds */
    public static final double MAX_JOB_TIME = 300.0;

    /** the xyz file that tinker will read from */
    public final TinkerXYZInputFile tinkerXYZInputFile;

    /** the key file that tinker will read from */
    public final TinkerKeyFile tinkerKeyFile;

    /**
     * Creates a job for minimizing some molecule.
     * @param molecule the molecule whose geometry is to be minimized
     * @param forcefield which forcefield to minimize on
     * @param extraKeywords any extra keywords that are desired (newlines required)
     */
    public TinkerMinimizationJob(Molecule molecule, Forcefield forcefield, String extraKeywords)
    {
        if ( molecule == null )
            throw new NullPointerException("molecule cannot be null when constructing a tinker minimization job");
        String keywords = "";
        if ( forcefield == Forcefield.AMOEBA )
            {
                keywords = DEFAULT_AMOEBA_KEYWORDS + extraKeywords;
       	        tinkerXYZInputFile = new TinkerXYZInputFile(molecule, Forcefield.AMOEBA);
            }
        else if ( forcefield == Forcefield.OPLS )
            {
                keywords = DEFAULT_OPLS_KEYWORDS + extraKeywords;
                tinkerXYZInputFile = new TinkerXYZInputFile(molecule, Forcefield.OPLS);
            }
       else
            throw new IllegalArgumentException("unrecognized forcefield type");
	    tinkerKeyFile = new TinkerKeyFile(keywords);
    }
    
    /** will run with standard keywords only */
    public TinkerMinimizationJob(Molecule molecule, Forcefield forcefield)
    {
        this(molecule,forcefield,"");
    }

    /**
     * Auto-selects a filename and runs the minimization calculation.
     * It is assumed that the TINKER binaries are in the path.
     * @return the result of the calculation
     */
    public TinkerMinimizationResult call()
    {
        // choose an appropriate base filename
        // try up to TINKER_MINIMIZATION_MAX_FILENAMES times to make a unique set of filenames
        String baseFilename = "";
        counting:
        for (int i=0; i < Settings.TINKER_MINIMIZATION_MAX_FILENAMES; i++)
            {
                // get a new ID number for this job
                int currentIndex = index.getAndIncrement();
                baseFilename = String.format("%s_tinker_minimization_job_%010d", Settings.HOSTNAME, currentIndex);

                // don't allow this choice of filenames if any files with this prefix already exist
                for ( File f : new File(Settings.TINKER_MINIMIZATION_JOB_DIRECTORY).listFiles() )
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
        if ( index.get() > Settings.TINKER_MINIMIZATION_MAX_FILENAMES )
            index.getAndSet(0);

	    // write input files to disk
        tinkerXYZInputFile.write( Settings.TINKER_MINIMIZATION_JOB_DIRECTORY + baseFilename + ".xyz" );
	    tinkerKeyFile.write     ( Settings.TINKER_MINIMIZATION_JOB_DIRECTORY + baseFilename + ".key" );

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
                        String runString = Settings.TINKER_MINIMIZATION_JOB_DIRECTORY + "run_tinker_minimization.sh " +
                                           Settings.TINKER_MINIMIZATION_JOB_DIRECTORY + " " + baseFilename;
                        //System.out.println(runString);
                        Process process = Runtime.getRuntime().exec(runString);
                        
                        while ( true )
                            {
                                // check if the process has terminated
                                GeneralThreadService.wait(500);
                               
                                /*try
                                    {
                                        TinkerMinimizationLogFile tempOutput = new TinkerMinimizationLogFile(Settings.TINKER_MINIMIZATION_JOB_DIRECTORY + baseFilename + ".out");
                                        if ( tempOutput.stringRepresentation.indexOf("Incomplete Convergence due to BadIntpln") > -1 )
                                            badGeometry = true;
                                    }
                                catch (Exception e) {}*/

                                boolean done = true;
                                try { exitValue = process.exitValue(); }
                                catch (Exception e) { done = false; }
                                //if ( done || badGeometry )
                                if ( done )
                                    break;

                                long now = System.currentTimeMillis();
                                elapsedTime = (now - startTime)/1000.0;
                                //if ( elapsedTime > 45.0 || badGeometry )
                                if ( elapsedTime > MAX_JOB_TIME )
                                    {
                                        process.destroy();
                                        break;
                                    }
                            }
                        long endTime = System.currentTimeMillis();
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
                        File[] files = new File(Settings.TINKER_MINIMIZATION_JOB_DIRECTORY).listFiles();
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
                        TinkerMinimizationLogFile errorOutput = new TinkerMinimizationLogFile(Settings.TINKER_MINIMIZATION_JOB_DIRECTORY + baseFilename + ".out");
                        String[] lines = errorOutput.stringRepresentation.split("\n");
                        int length = lines.length;
                        for (int i=Math.max(length-10,0); i < length; i++)
                            tail += lines[i] + "\n";
                    }
                catch (Exception e)
                    {
                    }
                throw new IllegalArgumentException("error code " + exitValue + " while runing tinker minimization job: " + baseFilename + "!\n" + tail);
            }

        // retrieve output XYZ File
	    TinkerXYZOutputFile xyzOutput = new TinkerXYZOutputFile(Settings.TINKER_MINIMIZATION_JOB_DIRECTORY + baseFilename + "_minimized.xyz");
        
        // retrieve file that contains energy
        TinkerMinimizationLogFile output = new TinkerMinimizationLogFile(Settings.TINKER_MINIMIZATION_JOB_DIRECTORY + baseFilename + ".out");
        
        // delete files after normal termination
        try
            {
                File[] files = new File(Settings.TINKER_MINIMIZATION_JOB_DIRECTORY).listFiles();
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
        return new TinkerMinimizationResult(output, elapsedTime, xyzOutput);
    }

    /** A class representing output files from a minimization call. Contains the molecule and log files produced from minimize call */
    public static class TinkerMinimizationResult implements Result, Serializable
    {
        /** for serialization */
        public static final long serialVersionUID = 1L;

	    /** The Tinker Minimiization Log File that is produced by minimize containing the energy of the optimized molecule */
	    public final TinkerMinimizationLogFile tinkerMinimizationLogFile;

	    /** The run time of the call in seconds. */
	    public final double elapsedTime;

	    /** The xyz file that is produced by a minimize call */
	    public final TinkerXYZOutputFile tinkerXYZOutputFile;

	    public TinkerMinimizationResult(TinkerMinimizationLogFile tinkerMinimizationLogFile, double elapsedTime, TinkerXYZOutputFile tinkerXYZOutputFile)
        {
	        this.tinkerMinimizationLogFile = tinkerMinimizationLogFile;
	        this.elapsedTime = elapsedTime;
	        this.tinkerXYZOutputFile = tinkerXYZOutputFile;
	    }
	
	    public Molecule getMolecule()
        {
	        return tinkerXYZOutputFile.molecule;
	    }

	    public double getEnergy()
        {
	        return tinkerMinimizationLogFile.energy;
	    }
	
	    @Override
	    public int hashCode()
        {
            return Objects.hash(tinkerMinimizationLogFile, elapsedTime, tinkerXYZOutputFile);
        }

        @Override
        public boolean equals(Object obj)
        {
            if ( obj == null )
                return false;
            if ( obj == this )
                return true;
            if ( !(obj instanceof TinkerMinimizationResult) )
                return false;

            TinkerMinimizationResult result = (TinkerMinimizationResult)obj;
            if ( Objects.equals(this.tinkerMinimizationLogFile, result.tinkerMinimizationLogFile) &&
                                this.elapsedTime == result.elapsedTime &&
                 Objects.equals(tinkerXYZOutputFile, result.tinkerXYZOutputFile) )
                return true;
            return false;
        }

        @Override
        public String toString()
        {
            return tinkerMinimizationLogFile.toString() + String.format(" time = %.3f s\n\n", elapsedTime) + tinkerXYZOutputFile.toString();
        }
    }
    
    /** A class representing the output from a minimization TINKER call */
    public static class TinkerMinimizationLogFile extends OutputFileFormat implements FileFormat
    {
	    /** The minimized energy in kcal/mol. */
	    public final Double energy;

        /** The final RMS gradient */
        public final Double gradient;

        /** The number of iterations */
        public final Integer iterations;

        /** the filename */
        public final String filename;

	    public TinkerMinimizationLogFile(String filename)
        {
            super(filename);
            this.filename = filename;

            // determine minimized energy
            Double readEnergy = null;
            Double readGradient = null;
            Integer readIterations = null;
            /*for ( int i=0; i < stringRepresentation.length(); i++ )
                {
                    String currentLetter = stringRepresentation.substring(i,i+1);
                    char thisChar = currentLetter.charAt(0);
                    int code = (int)thisChar;
                    System.out.println(currentLetter + " (" + code + ")");
                }*/
            for (String line : stringRepresentation.split("\n"))
                {
                    //System.out.println(index + " " + line);
                    String[] fields = line.trim().split("\\s+");
                    if ( line.indexOf("Final Function Value") > -1 )
                        readEnergy = Double.valueOf(fields[4]);
                    else if ( line.indexOf("Final RMS Gradient") > -1 )
                        {
                            try { readGradient = Double.valueOf(fields[4]); }
                            catch (Exception e) { readGradient = 999.9; }
                            //System.out.println(readGradient + " : " + line);
                        }
                    else 
                        {
                            try
                                {
                                    readIterations = Integer.valueOf(fields[0]);
                                }
                            catch (NumberFormatException e)
                                {
                                }
                        }
                }

            if ( readEnergy == null )
                {
                    String exceptionString = "Energy not found!\n";
                    int start = Math.max(0, fileContents.size()-5);
                    for (int i=start; i < fileContents.size(); i++)
                        exceptionString += fileContents.get(i) + "\n";
                    throw new IllegalArgumentException(exceptionString);
                }
            if ( readGradient == null )
                throw new IllegalArgumentException("Gradient not found!");
            if ( readIterations == null )
                throw new IllegalArgumentException("Number of iterations not found!");
            this.energy = readEnergy;
            this.gradient = readGradient;
            this.iterations = readIterations;
        }

	    @Override
	    public int hashCode()
        {
            return Objects.hash(stringRepresentation, energy);
        }

        @Override
        public boolean equals(Object obj)
        {
            if ( obj == null )
                return false;
            if ( obj == this )
                return true;
            if ( !(obj instanceof TinkerMinimizationLogFile) )
                return false;

            TinkerMinimizationLogFile file = (TinkerMinimizationLogFile)obj;
            if ( Objects.equals(this.stringRepresentation, file.stringRepresentation) &&
                 Objects.equals(this.energy, file.energy ) )
                return true;
            return false;
        }

        @Override
        public String toString()
        {
            return String.format("TinkerMinimizationLogFile: energy = %.3f", energy);
        }
    }

    @Override
    public int hashCode()
    {
        return Objects.hash(tinkerXYZInputFile, tinkerKeyFile);
    }

    @Override
    public boolean equals(Object obj)
    {
        if ( obj == null )
            return false;
        if ( obj == this )
            return true;
        if ( !(obj instanceof TinkerMinimizationJob) )
            return false;

        TinkerMinimizationJob job = (TinkerMinimizationJob)obj;
        if ( Objects.equals(this.tinkerXYZInputFile, job.tinkerXYZInputFile) &&
             Objects.equals(this.tinkerKeyFile,      job.tinkerKeyFile) )
            return true;
        return false;
    }

    @Override
    public String toString()
    {
        return "TinkerJob:\n" + tinkerXYZInputFile.stringRepresentation + "\n" + tinkerKeyFile.stringRepresentation;
    }

    /** For testing. */
    public static void main(String[] args)
    {
        DatabaseLoader.go();
        List<ProtoAminoAcid> sequence = ProtoAminoAcidDatabase.getSpecificSequence("arg","met","standard_ala","gly","d_proline", "gly", "phe", "val", "hd", "l_pro");
        Peptide peptide = PeptideFactory.createPeptide(sequence);
        
        for (int i=0; i < peptide.sequence.size(); i++)
            {
                peptide = BackboneMutator.mutateOmega(peptide, i);
                peptide = BackboneMutator.mutatePhiPsi(peptide, i);
                peptide = RotamerMutator.mutateChis(peptide, i);
            }
        peptide = PeptideFactory.setHairpinAngles(peptide);

        TinkerMinimizationJob job = new TinkerMinimizationJob(peptide, Forcefield.AMOEBA, "maxiter 2000\n\n");
        TinkerMinimizationJob.TinkerMinimizationResult result = job.call();
        Molecule minimizedMolecule = result.tinkerXYZOutputFile.molecule;
        Peptide newPeptide = peptide.setPositions(minimizedMolecule);

        double time = result.elapsedTime;
        TinkerMinimizationJob.TinkerMinimizationLogFile logFile = result.tinkerMinimizationLogFile;
        double energy = logFile.energy;
        double gradient = logFile.gradient;
        int iterations = logFile.iterations;

        for (int i=0; i < peptide.contents.size(); i++)
            System.out.println(peptide.contents.get(i).toFullString() + "   " + newPeptide.contents.get(i).toFullString());
        System.out.printf("\ntime %.1f s   energy %.2f kcal   RMS gradient %.4f   %d iterations\n", time, energy, gradient, iterations);
    }
}
