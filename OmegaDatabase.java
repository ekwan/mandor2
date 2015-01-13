import java.util.*;
import java.io.*;
import org.apache.commons.math3.distribution.NormalDistribution;
import com.google.common.collect.*;

/**
 * This singleton stores the backbone-dependent omega data.  Omega is defined as the
 * the dihedral angle of the amide bond.  Random numbers are drawn from a normal
 * distribution centered on the given mean and standard deviation in the database.
 * A new NormalDistribution and associated apache RandomGenerator are generated
 * for each draw.  To be more multithreading efficient, we could delegate the random number
 * generator to the calling class, but we're leaving it this way for now.
 *
 * The library tells you how the omega(+1) dihedral depends on the current
 * torsion, psi(0), and the following torsion, phi(+1).  Here is the diagram
 * from the library:<p>
 * <p><pre>
 *          omega(0)    phi(0)     psi(0)  omega(+1)   phi(+1)    psi(+1) omega(+2)<p>
 *              -  N( 0)  -  Ca( 0)  -  C( 0)  -  N(+1)  -  Ca(+1)  -  C(+1)  -    <p>
 *                       current residue               following residue           <p>
 * <p></pre>
 * The library contains some global data starting with "All" but that data are ignored
 * because we assume we're given the amino acids on either side of the torsion.  There's
 * a confusing convention that omega precedes an amino acid in the N to C direction, but
 * we follow it anyways.
 */
public class OmegaDatabase implements Singleton
{
    /** Not instantiable. */
    public OmegaDatabase()
    {
        throw new IllegalArgumentException("not instantiable");
    }
    
    /**
     * Stores all of the omega data.
     * The String key represents the type of data (ResTypeGroup: e.g., All, All_nonxpro, etc.)
     * The keys: first pair is phi(+1),psi(0), second pair is the given mean and standard deviation
     */
    private static final Map<String,Map<Pair<Double,Double>,Pair<Double,Double>>> DATABASE;

    /** Static initializer. */
    static
    {
        HashMap<String,Map<Pair<Double,Double>,Pair<Double,Double>>> database = new HashMap<>();
        String filename = Settings.OMEGA_DATA_FILENAME;

        // read data from file
        Scanner thisFile = null;
        try
            {
                thisFile = new Scanner(new FileReader(filename));
                String lastResidue = null;
                Map<Pair<Double,Double>,Pair<Double,Double>> currentData = new LinkedHashMap<>();
                while (thisFile.hasNextLine())
                    {
                        String currentLine = thisFile.nextLine().trim();

                        // ignore header lines and "all-residue" lines
                        if ( currentLine.startsWith("@") || currentLine.startsWith("#") ||
                             currentLine.startsWith("All") )
                            continue;
                        
                        // parse fields
                        String[] fields = currentLine.split("\\s+");
                        String currentResidue = fields[0];
                        if ( lastResidue == null )
                            lastResidue = currentResidue;
                        
                        // if this is a new residue type, load the stored data into the database
                        // assumes file ends with a line of data
                        if ( ! lastResidue.equals(currentResidue) )
                            {
                                database.put(lastResidue,currentData);
                                currentData = new LinkedHashMap<>();
                            }
                        lastResidue = currentResidue;

                        // record the current data
                        Double phi = Double.valueOf(fields[1]);
                        Double psi = Double.valueOf(fields[2]);
                        Pair<Double,Double> anglePair = new Pair<>(phi,psi);
                        Double mean = Double.valueOf(fields[5]);
                        Double stdev = Double.valueOf(fields[6]);
                        Pair<Double,Double> distPair = new Pair<>(mean,stdev);
                        currentData.put(anglePair,distPair);
                    }

                // load last entry into database
                database.put(lastResidue, currentData);
            }
        catch (IOException e)
            {
                e.printStackTrace();
                System.exit(1);
            }
        finally
            {
                if ( thisFile != null )
                    thisFile.close();
            }
        DATABASE = ImmutableMap.copyOf(database);
    }

    /**
     * Returns a normal distribution for the omega torsion of residue2, based on the
     * identities of residue1 and residue2 and the enclosing phi angles, psi0 and
     * phi1.
     * @param aa1 the amino acid before the omega torsion (N to C direction)
     * @param psi0 the psi of the preceding amino acid aa1 in degrees
     * @param aa2 the amino acid after the omega torsion (N to C direction)
     * @param phi1 the phi of the following amino acid aa2 in degrees
     * @return the normal distribution for the residue2 omega torsion
     */
    public static NormalDistribution getNormalDistribution(AminoAcid aa1, double psi0, AminoAcid aa2, double phi1)
    {
        String key = determineKey(aa1, aa2);
        Map<Pair<Double,Double>,Pair<Double,Double>> thisEntry = DATABASE.get(key);
        if ( thisEntry == null )
            throw new NullPointerException("should have found something in the database!");
        
        // data are stored as phi1, psi0:
        // @ ResTypeGroup  Phi(+1) Psi(0)  S Num     mW(+1) sW(+1)
        Pair<Double,Double> key2 = snapToGrid(new Pair<>(phi1,psi0));
        //System.out.println(key2);
        //for (Pair<Double,Double> p : thisEntry.keySet() )
        //    System.out.println( p.toString() + " : " + (p.getFirst() - key2.getFirst()) + " , " + (p.getSecond() - key2.getSecond()) );
        Pair<Double,Double> dist = thisEntry.get(key2);
        if ( dist == null )
            {
                String errorString = String.format("unexpected error in getNormalDistribution: residue1=%s residue2=%s", 
                                                   aa1.shortName, aa2.shortName);
                throw new NullPointerException(errorString);
            }
        return new NormalDistribution(dist.getFirst(), dist.getSecond());
    }

    /**
     * Returns a random omega given psi0 of the preceding residue1 and the phi1 of the
     * following residue2.  Nulls are not allowed.
     * @param aa1 the amino acid before the omega torsion (N to C direction)
     * @param psi0 the psi of the preceding amino acid aa1 in degrees
     * @param aa2 the amino acid after the omega torsion (N to C direction)
     * @param phi1 the phi of the following amino acid aa2 in degrees
     * @return a random, normally-distributed value for the omega of residue2 based on the database values
     */
    public static double getOmega(AminoAcid aa1, double psi0, AminoAcid aa2, double phi1)
    {
        if ( aa1 == null || aa2 == null )
            throw new IllegalArgumentException("nulls are not allowed");
        if ( psi0 < -180.0 || psi0 > 180.0 )
            throw new IllegalArgumentException("psi0 out of range");
        if ( phi1 < -180.0 || phi1 > 180.0 )
            throw new IllegalArgumentException("phi1 out of range");
        NormalDistribution dist = getNormalDistribution(aa1, psi0, aa2, phi1);
        return dist.sample();
    }

    /**
     * Returns the database key for this combination of residues.
     * Amino acid 1 precedes amino acid 2 in the N to C direction.
     * @param aa1 the residue before the omega1 torsion
     * @param aa2 the residue after the omega1 torsion
     * @return the database key
     */
    private static String determineKey(AminoAcid aa1, AminoAcid aa2)
    {
        String key = "";

        if ( aa1 == AminoAcid.ILE || aa1 == AminoAcid.VAL )
            key = "IleVal";
        else if ( aa1 == AminoAcid.GLY )
            key = "Gly";
        else if ( aa1.isProline() )
            key = "Pro";
        else
            key = "NonPGIV";

        if ( aa2.isProline() )
            {
                // the torsion precedes a proline
                key = key + "_xpro";
            }
        else
            {
                // the torsion does not precede a proline
                key = key + "_nonxpro";
            }
        //System.out.println("key is " + key);
        return key;
    }

    /**
     * rounds phi and psi to nearest multiple of 10
     * @param inputPair e.g. phi1,psi0 (-179.0, -178.0) 
     * @return the nearest pair on the grid e.g., phi, psi (-180.0, 180.0)
     */
    private static Pair<Double,Double> snapToGrid(Pair<Double,Double> inputPair)
    {
        double phi = inputPair.getFirst();
        double psi = inputPair.getSecond();
        if (phi > 180.0 || phi < -180.0 || psi > 180.0 || psi < -180.0)
            throw new IllegalArgumentException("psi and phi must be between -180 and 180");

        double phi_rounded = NonRotamericLibrary.roundTo10(phi);
        double psi_rounded = NonRotamericLibrary.roundTo10(psi);
        
        if (phi_rounded == 180.0)
            phi_rounded = -180.0;

        if (psi_rounded == 180.0)
            psi_rounded = -180.0;

        // return the appropriate data
        return new Pair<Double,Double>(phi_rounded,psi_rounded);
    }

    /** Forces the static initializer to run. */
    public static void load()
    {
        // counts how many lines of data were read
        int entries = 0;
        for (String s : DATABASE.keySet())
            {
                for (Pair<Double,Double> p : DATABASE.get(s).keySet())
                    entries++;
            }
        System.out.printf("Omega dihedral angle database loaded with %d entries.\n", entries);
    }

    /** Tester class. */
    public static void main(String[] args)
    {
        AminoAcid aa1 = AminoAcid.ARG;
        AminoAcid aa2 = AminoAcid.PHE;
        double psi0 = 120.0;
        double phi1 = -120.0;
        NormalDistribution dist = getNormalDistribution(aa1, psi0, aa2, phi1);
        System.out.println(dist.getMean());
        System.out.println(dist.getStandardDeviation());
        System.out.println(getOmega(aa1, psi0, aa2, phi1));
    }
}
