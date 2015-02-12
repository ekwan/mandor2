import java.util.*;
import java.io.*;
import com.google.common.collect.*;

/**
 * Holds variables that control the behavior of the program.
 */
public class Settings implements Immutable, Singleton
{
    // ProtoAminoAcidDatabase Settings
    
        /** the directory containing the templates for all the amino acids */
        public static final String PROTOAMINOACID_DIRECTORY;

    // Rotamer Library Database Settings

        /** the directory containing the Dunbrack backbone-dependent rotamer data */
        public static final String ROTAMER_LIBRARY_DIRECTORY = "rotamer_library/";

        /** the probability below which rotamers will be ignored */
        public static final double ROTAMER_LIBRARY_THRESHOLD = 0.01;

    // Ramachandran Database Settings

        // the directory where the backbone-dependent Ramachandran data are stoed
        public static final String RAMACHANDRAN_DIRECTORY = "databases/";
        
        // the prefix of all the backbone-dependent Ramachandran data
        public static final String RAMACHANDRAN_DATA_PREFIX = "NDRD_split";

    // Omega Database Settings

        /** the file containing the omega data */
        public static final String OMEGA_DATA_FILENAME = "databases/omegaCDL_OmegaBetweenAsPhi1Psi0_KernRegr_v1.3.1_Aug12-2011.txt";

    // OPLScalculator Settings

        /** the file containing the OPLS forcefield file */
        public static final String OPLS_FORCEFIELD_FILENAME = "amino_acids/oplsaal.prm";

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

    // Tinker General Parameters

        /** The directory where the tinker job setup script is. */
        public static final String TINKER_SETUP_DIRECTORY;

        /** The directory where the tinker job setup script should make some directories. */
        public static final String TINKER_TARGET_DIRECTORY;

        /** The directory where the tinker executables are. */
        public static final String TINKER_PROGRAM_DIRECTORY;

    // Tinker Minimization Job Parameters

        /** where tinker minimization jobs will be run */
        public static final String TINKER_MINIMIZATION_JOB_DIRECTORY;
    
        /** how many temporary filenames are available for running Tinker minimization jobs */
        public static final int TINKER_MINIMIZATION_MAX_FILENAMES = 20000;

    // Tinker Analysis Job Parameters
    
        /** where tinker analysis jobs will be run */
        public static final String TINKER_ANALYSIS_JOB_DIRECTORY;
        
        /** how many temporary filenames are available for running Tinker minimization jobs */
        public static final int TINKER_ANALYSIS_MAX_FILENAMES = 20000;

    // Steric Energy Settings

        /** cutoff distance in angstroms for calculating steric energies */
        public static final double CUTOFF_DISTANCE = 6.0;

    // Geometry Settings

        /** atoms farther this apart in angstroms are considered not hydrogen bonded */
        public static final double MAXIMUM_HBOND_DISTANCE = 2.20;

        /** minimum H-bond angle in degrees */
        public static final double MINIMUM_HBOND_ANGLE = 120.0;

        /** polar elements in a hydrogen bond */
        public static final Set<Element> HBOND_ELEMENTS = ImmutableSet.of(Element.OXYGEN, Element.NITROGEN, Element.SULFUR);

        /** atoms closer than this distance in angstroms are considered to be clashing */
        public static final double MINIMUM_INTERATOMIC_DISTANCE = 0.90;

    /** static initializer */
    static
    {
        String temp = "";

        // set hostname
        try { temp = java.net.InetAddress.getLocalHost().getHostName(); } catch (Exception e) {}
        if ( temp.length() == 0 )
            temp = "localhost";
        FULL_HOSTNAME = temp;

        if ( FULL_HOSTNAME.length() > 0 )
            temp = FULL_HOSTNAME.split("\\.")[0];
        else
            temp = "";
        HOSTNAME = temp;

        // set number of threads
        int tempThreads = Runtime.getRuntime().availableProcessors();
        if (HOSTNAME.startsWith("enj"))
            tempThreads = 12;
        else if ( HOSTNAME.startsWith("dae"))
            tempThreads = 8;
        else if ( HOSTNAME.startsWith("holy"))
            tempThreads = 64;
        else if ( HOSTNAME.equals("localhost") )
            tempThreads = 2;

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

        // set up ProtoAminoAcid directory
        PROTOAMINOACID_DIRECTORY = WORKING_DIRECTORY + "amino_acids/";

        // set up tinker directories
        TINKER_SETUP_DIRECTORY = WORKING_DIRECTORY + "tinker_jobs/";
        if ( HOSTNAME.startsWith("dae") || HOSTNAME.startsWith("enj") || HOSTNAME.startsWith("holy") )
            {
                TINKER_TARGET_DIRECTORY = "/dev/shm/";
                TINKER_PROGRAM_DIRECTORY = "/n/jacobsen_lab/ekwan/tinker7/bin/";
            }
        else if ( HOSTNAME.equals("localhost") )
            {
                TINKER_TARGET_DIRECTORY = WORKING_DIRECTORY + "tinker_jobs/";
                TINKER_PROGRAM_DIRECTORY = "/Users/ekwan/tinker/bin/";
            }
        else
            throw new IllegalArgumentException("recognized computer");
        String runString = Settings.TINKER_SETUP_DIRECTORY + "tinker_setup.sh " + TINKER_TARGET_DIRECTORY + " "
                           + PROTOAMINOACID_DIRECTORY + " " + TINKER_SETUP_DIRECTORY + " " + TINKER_PROGRAM_DIRECTORY;
        try
            {
                Process process = Runtime.getRuntime().exec(runString);
                process.waitFor();
                if ( process.exitValue() != 0 )
                    throw new IllegalArgumentException("tinker setup terminated abormally, exit code: " + process.exitValue());
            }
        catch (Exception e)
            {
                e.printStackTrace();
                System.exit(1);
            }

        // for TinkerMinimizationJobs
        TINKER_MINIMIZATION_JOB_DIRECTORY = TINKER_TARGET_DIRECTORY + "tinker_minimization_jobs/";

	    // for TinkerAnalysisJobs
	    //temp = WORKING_DIRECTORY + "tinker_analysis_jobs/";
        //temp = "/dev/shm/tinker_analysis_jobs/";
	    //if ( PLATFORM == Platform.DOS )
        //    temp = temp.replace("/","\\");
        TINKER_ANALYSIS_JOB_DIRECTORY = TINKER_TARGET_DIRECTORY + "tinker_analysis_jobs/";

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
        throw new IllegalArgumentException("not instantiable!");
    }

    /** For testing. */
    public static void main(String[] args)
    {
        System.out.println("hello from Settings");
    }
}
