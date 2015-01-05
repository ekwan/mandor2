import java.util.*;
import java.io.*;

/**
 * Holds variables that control the behavior of the program.
 */
public class Settings implements Immutable, Singleton
{
    // ProtoAminoAcidLibrary Settings
    
        /** the directory containing the templates for all the amino acids */
        public static final String PROTO_AMINO_ACID_LIBRARY_DIRECTORY = "input_files/";

    // Rotamer Library Settings

        /** the directory containing the Dunbrack backbone-dependent rotamer data */
        public static final String ROTAMER_LIBRARY_DIRECTORY = "rotamer_library/";

        /** the probability below which rotamers will be ignored */
        public static final double ROTAMER_LIBRARY_THRESHOLD = 0.01;

    // Ramachandran Library Settings

        /** the file containing the backbone-dependent Ramachandran data */
        public static final String RAMACHANDRAN_DATA_FILENAME = "ramachandran/NDRD_Tonly.txt.gz";

    // Omega Library Settings

        /** the file containing the omega data */
        public static final String OMEGA_DATA_FILENAME = "omega/omegaCDL_OmegaBetweenAsPhi1Psi0_KernRegr_v1.3.1_Aug12-2011.txt";

    // General Settings

        /** hostname (e.g., enj02.rc.fas.harvard.edu) */
        public static final String FULL_HOSTNAME;

        /** first field of hostname */
        public static final String HOSTNAME;

        /** number of available threads */
        public static final int NUMBER_OF_THREADS;

        /** type of system */
        public static final Platform PLATFORM;
        
        /** current working directory */
        public static final String WORKING_DIRECTORY;

        /** the main class name */
        public static final String MAIN_CLASS;

        public enum Platform
        {
            DOS, LINUX;
        }

    // Omnisol Job Parameters

        /** location of Omnisol program */
        public static final String OMNISOL_LOCATION;

        /** input header for Omnisol chloroform calculation */
        public static final String OMNISOL_HEADER_CHLOROFORM;

        /** input header for Omnisol water calculation */
        public static final String OMNISOL_HEADER_WATER;

        /** how many temporary filenames are available for running Omnisol jobs */
        public static final int OMNISOL_MAX_FILENAMES = 10000;

    // Tinker Minimization Job Parameters

        /** where tinker minimization jobs will be run */
        public static final String TINKER_MINIMIZATION_JOB_DIRECTORY;
    
        /** standard keywords to use for every tinker minimization job */
        public static final String TINKER_MINIMIZATION_STANDARD_KEYWORDS = "parameters amoebapro13.prm\nwriteout 200\n\n";

        /** how many temporary filenames are available for running Tinker minimization jobs */
        public static final int TINKER_MINIMIZATION_MAX_FILENAMES = 20000;

        /** object to synchronize on for copying files into /dev/shm */
        public static final Object SHM_LOCK = new Object();

    // Tinker Analysis Job Parameters
    
        /** where tinker analysis jobs will be run */
        public static final String TINKER_ANALYSIS_JOB_DIRECTORY;
        
        /** where gaussian analysis jobs will be run */
        public static final String GAUSSIAN_ANALYSIS_JOB_DIRECTORY;

        /** standard keywords to use for every tinker analysis job */
        public static final String TINKER_ANALYSIS_STANDARD_KEYWORDS = "parameters amoebapro13.prm\n\n";

        /** how many temporary filenames are available for running Tinker minimization jobs */
        public static final int TINKER_ANALYSIS_MAX_FILENAMES = 20000;

    // Forcefield Job Parameters

        /** where forcefield jobs will be run */
        public static final String FORCEFIELD_JOB_DIRECTORY;

        /** how many temporary filenames are available for running forcefield jobs */
        public static final int FORCEFIELD_JOB_MAX_FILENAMES = 10000;

    // ProtoAminoAcid settings

        /** location of amino acid template files */
        public static final String AMINO_ACID_DIRECTORY = "amino_acids/";

    // Steric Energy Settings

        /**
         * the threshold distance in angstroms for an interatomic distance to be
         * considered too small (ignores directly connected atoms)
         */
        public static final double MINIMUM_DISTANCE = 1.00;

        /** cutoff distance in angstroms for calculating steric energies */
        public static final double CUTOFF_DISTANCE = 6.0;

    /** static initializer */
    static
    {
        String temp = "";

        // set hostname
        try { temp = java.net.InetAddress.getLocalHost().getHostName(); } catch (Exception e) {}
        FULL_HOSTNAME = temp;

        if ( FULL_HOSTNAME.length() > 0 )
            temp = FULL_HOSTNAME.split("\\.")[0];
        else
            temp = "";
        HOSTNAME = temp;

        // set number of threads
        int tempThreads = Runtime.getRuntime().availableProcessors();
        if (HOSTNAME.startsWith("enj"))
            tempThreads = 13;
        else if ( HOSTNAME.startsWith("dae"))
            tempThreads = 9;
        else if ( HOSTNAME.startsWith("holy"))
            tempThreads = 65;
        
        //NUMBER_OF_THREADS=12;
        NUMBER_OF_THREADS = tempThreads;

        // detect platform
        temp = System.getProperty("os.name").toLowerCase();
        if ( temp.indexOf("win") >= 0 )
            PLATFORM = Platform.DOS;
        else if ( temp.indexOf("nix") >= 0 || temp.indexOf("nux") >= 0 || temp.indexOf("aix") >= 0 || temp.indexOf("mac") >= 0 )
            PLATFORM = Platform.LINUX;
        else
            throw new IllegalArgumentException("Unsupported operating system: " + temp);

        // get working directory
        temp = System.getProperty("user.dir") + "/";
        if ( PLATFORM == Platform.DOS )
            temp = temp.replace("/","\\");
        WORKING_DIRECTORY = temp;

        // set OMNISOL_HEADER_CHLOROFORM
        temp =        "SM5.0R\n";
        temp = temp + "& IOFR=1.4459 ALPHA=0.15 BETA=0.02 GAMMA=38.39\n";
        temp = temp + "& FACARB=0.00 FEHALO=0.75 SOLVNT=GENORG\n";
        temp = temp + "peptide (solvent : chloroform)\n\n";
        OMNISOL_HEADER_CHLOROFORM = temp;

        //set OMNISOL_HEADER_WATER
        temp = "SM5.0R SOLVNT=WATER\n";
        temp = temp + "peptide (aqueous)\n\n";
        OMNISOL_HEADER_WATER = temp; 

        // set OMNISOL_LOCATION
        temp = WORKING_DIRECTORY + "omnisol/";
        if ( PLATFORM == Platform.DOS )
            temp = temp.replace("/","\\");
        OMNISOL_LOCATION = temp;

        // for TinkerMinimizationJobs
        //temp = WORKING_DIRECTORY + "tinker_minimization_jobs/";
        temp = "/dev/shm/tinker_minimization_jobs/";
        if ( PLATFORM == Platform.DOS )
            temp = temp.replace("/","\\");
        TINKER_MINIMIZATION_JOB_DIRECTORY = temp;

	    // for TinkerAnalysisJobs
	    //temp = WORKING_DIRECTORY + "tinker_analysis_jobs/";
        temp = "/dev/shm/tinker_analysis_jobs/";
	    if ( PLATFORM == Platform.DOS )
            temp = temp.replace("/","\\");
        TINKER_ANALYSIS_JOB_DIRECTORY = temp;

        // for GaussianAnalysisJobs
        temp = WORKING_DIRECTORY + "gaussian_analysis_jobs/";
        GAUSSIAN_ANALYSIS_JOB_DIRECTORY = temp;

        // for ForcefieldJobs
        temp = WORKING_DIRECTORY + "forcefield_jobs/";
        if ( PLATFORM == Platform.DOS )
            temp = temp.replace("/","\\");
        FORCEFIELD_JOB_DIRECTORY = temp;

        // print out some information
        //System.out.println(String.format("Mandor hostname is %s (%d cores available).", HOSTNAME, NUMBER_OF_THREADS));
    
        // set the main class name
        StackTraceElement[] stack = Thread.currentThread().getStackTrace();
        StackTraceElement main = stack[stack.length - 1];
        MAIN_CLASS = main.getClassName();
    }

    /** not instantiable */
    private Settings()
    {
    }

    public static void main(String[] args)
    {
        System.out.println("test");
    }
}
