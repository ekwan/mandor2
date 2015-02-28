import java.util.*;
import java.util.concurrent.*;
import java.util.concurrent.atomic.*;
import com.google.common.collect.*;
import org.apache.commons.io.*;
import java.io.*;

/**
 * This class collects together some methods for doing Monte Carlo minimizations of peptide structures.
 */
public abstract class MonteCarloJob implements WorkUnit
{
    /** The starting structure. This should be minimized and have an energy breakdown. */
    public final Peptide startingPeptide;
    
    /** The amount to decrease the temperature by on each iteration. */
    public final double deltaAlpha;

    /** How many iterations to perform. */
    public final int maxIterations;

    public MonteCarloJob()
    {
        startingPeptide = null;
        deltaAlpha = 0.0;
        maxIterations = 0;
    }

    /** Constructs a generic MonteCarloJob. */
    public MonteCarloJob(Peptide startingPeptide, double deltaAlpha, int maxIterations)
    {
        this.startingPeptide = startingPeptide;
        this.deltaAlpha = deltaAlpha;
        this.maxIterations = maxIterations;
    }

    /**
     * Makes a random mutation of some sort to peptide and returns the new Peptide.
     * @param peptide the starting structure
     * @return the final structure
     */
    public abstract Peptide mutate(Peptide peptide);

    /**
     * Decides whether we should accept the new structure based on the
     * modified Metropolis criterion as described by Holm and Sanders in
     * PROTEINS: Structure, Function, and Genetics, 1992, 14, 213-223.
     * This method is thread safe.
     * @param oldPeptide the starting structure
     * @param newPeptide the candidate structure
     * @param currentAlpha the current inverse temperature of the simulation (1/K)
     */
    public static boolean acceptChange(Peptide oldPeptide, Peptide newPeptide, double currentAlpha)
    {
        double oldEnergy = oldPeptide.energyBreakdown.totalEnergy;
        double newEnergy = newPeptide.energyBreakdown.totalEnergy;
        double deltaE    = newEnergy - oldEnergy;
        if ( deltaE < 0.0 )
            return true;
        else
            {
                double probability = 1.0 / ( 1.0 + Math.exp(currentAlpha * deltaE) );
                double draw = ThreadLocalRandom.current().nextDouble();
                if ( draw <= probability )
                    return true;
                return false;
            }
    }

    /**
     * Runs a generic Monte Carlo simulation.
     * @param maxSize the number of best results to keep
     * @return the list of the best structures
     */
    public List<Peptide> minimize(int maxSize)
    {
        // this is the current state of the Monte Carlo simulation
        Peptide currentPeptide = startingPeptide;
        Peptide lastPeptide = startingPeptide;

        // this is where we will store the best peptides
        PeptideList bestResults = new PeptideList(maxSize);
        boolean usingLast = true;
        bestResults.add(currentPeptide);
        
        // perform the Metropolis Monte Carlo
        double currentAlpha = 0.0;
        for (int i=0; i < maxIterations; i++)
            {
                // abort if kill file found
                if ( new File("kill.txt").isFile() )
                    {
                        System.out.println("Kill file found.  Aborting.");
                        break;
                    }

                // make a mutation
                Peptide candidate = mutate(currentPeptide);

                // minimize the candidate
                //candidate = minimizeSingle(candidate);

                // apply modified Metropolis criterion
                boolean isAccepted = acceptChange(currentPeptide, candidate, currentAlpha);
                if ( isAccepted )
                    {
                        lastPeptide = currentPeptide;
                        usingLast = false;
                        currentPeptide = candidate;
                        bestResults.add(currentPeptide);
                    }

                // write information about this iteration
                System.out.printf("[ %s ] Iter %d of %d (alpha = %.6f, E = %.2f, bestE = %.2f)\n", new Date().toString(),
                                  i+1, maxIterations, currentAlpha, candidate.energyBreakdown.totalEnergy, bestResults.getBestEnergy());

                // update alpha
                currentAlpha = currentAlpha + deltaAlpha;
            }
        
        return bestResults.getList();
    }

    /**
     * Appends a string to a file.
     * @param logFile the file to write to
     * @param reportString the string to be appended
     */
    public static void appendStringToFile(File logFile, String reportString)
    {
        try
            {
                FileUtils.writeStringToFile(logFile, reportString + "\n", "UTF-8", true);
            }
        catch (Exception e)
            {
                System.out.println("Error writing to log file " + logFile.getName() + "!");
                e.printStackTrace();
            }
    }

    /** Runs the minimization. */
    public abstract MonteCarloResult call();

    /**
     * Represents a list of peptides with a maximum size.  This object is not thread safe.
     */
    public static class PeptideList implements Serializable
    {
        /** for serialization */
        public static final long serialVersionUID = 1L;

        /** Maximum size of the peptide list. */
        private final int maxSize;

        /** The internal list of peptides. */
        private final LinkedList<Peptide> list;

        /**
         * Creates a list of peptides with a maximum size.  Only keeps the lowest energy results.
         * @param maxSize the maximum size of the list
         */
        public PeptideList(int maxSize)
        {
            this.maxSize = maxSize;
            list = new LinkedList<>();
        }

        /**
         * Attempts to add a peptide to the list.  If this peptide is low in energy it will be retained.
         * @param peptide the peptide to add
         */
        public void add(Peptide peptide)
        {
            if ( peptide == null )
                throw new NullPointerException("nulls not supported");
            if ( peptide.energyBreakdown == null || peptide.energyBreakdown == EnergyBreakdown.BLANK )
                throw new IllegalArgumentException("peptide EnergyBreakdown not set!");
            
            double thisEnergy = peptide.energyBreakdown.totalEnergy;
            
            if ( list.size() < maxSize )
                {
                    //System.out.println("list is not full, added: " + peptide.energyBreakdown.totalEnergy);
                    list.add(peptide);
                    Collections.sort(list);
                }
            else
                {
                    double worstEnergy = list.get(list.size()-1).energyBreakdown.totalEnergy;
                    if ( thisEnergy < worstEnergy )
                        {
                            list.removeLast();
                            list.add(peptide);
                            Collections.sort(list);
                            //for (Peptide p : list)
                            //    System.out.println(p.energyBreakdown.totalEnergy);
                        }
                }
            
            //if ( peptide == list.get(0) )
            //    System.out.printf("*** New Best: %.2f, dropped by %.2f ***", thisEnergy, thisEnergy - list.get(1).energyBreakdown.totalEnergy );
        }

        public List<Peptide> getList()
        {
            return ImmutableList.copyOf(list);
        }

        public double getBestEnergy()
        {
            if ( list == null || list.size() == 0 )
                return 0.0;
            return list.get(0).energyBreakdown.totalEnergy;
        }

        public Peptide get(int i)
        {
            return list.get(i);
        }
    }
}
