import java.util.*;
import java.util.zip.GZIPInputStream;
import java.io.*;
import com.google.common.collect.*;
import com.google.common.primitives.*;
import java.util.concurrent.*;

/**
 * This singleton stores data on how (phi,psi) depends on the neighboring residue identity based on the Dunbrack data.<p>
 * Reference: Daniel Ting, Guoli Wang, Maxim Shapovalov, Rajib Mitra, Michael I. Jordan, Roland L. Dunbrack, Jr.
 * <u>Neighbor-dependent Ramachandran probability distributions of amino acids developed from a</u>
 * <u>hierarchical Dirichlet process model.</u><i> PLOS Comp. Biol.</i><em> (April 2010)</em>.<p><p>
 *
 * The libraries are now read in parallel.  Data are not available for Xxx-CPR or CPR-Xxx in the main files, though
 * there are CPR-Xxx and Xxx-CPR data in the CPR file.  (CPR = cis-proline.)  I could probably figure out a way to
 * deal with that, but for now, I just treat it like it's trans-proline (TPR).
 *
 * Use this command to check memory usage:<p>
 * <code>jps | grep Rama | awk '{print $1}' | xargs jmap -histo:live | awk '{if ($1 == "Total" || NR < 20) {print}}'
 */
public class RamachandranDatabase implements Singleton
{
    /** The Ramachandran data. */
    private static final Map<CustomKey,PreDistribution> DATABASE;

    /** read all the Ramachandran data */
    static
    {
        // temporary copy of the database
        Map<CustomKey,PreDistribution> tempDatabase = new ConcurrentHashMap<>();

        // load entries in parallel
        String databaseDirectory = Settings.WORKING_DIRECTORY + Settings.RAMACHANDRAN_DIRECTORY;
        List<Future<Result>> futures = new ArrayList<>();
        for (File file : new File(databaseDirectory).listFiles())
            {
                String filename = Settings.WORKING_DIRECTORY + Settings.RAMACHANDRAN_DIRECTORY + file.getName();
                if ( filename.endsWith("gz") && file.getName().startsWith(Settings.RAMACHANDRAN_DATA_PREFIX) && filename.indexOf("header") == -1 )
                    {
                        RamachandranUnit unit = new RamachandranUnit(filename, tempDatabase);
                        Future<Result> f = GeneralThreadService.submit(unit);
                        futures.add(f);
                    }
            }
        GeneralThreadService.silentWaitForFutures(futures);

        // make immutable copy
        DATABASE = ImmutableMap.copyOf(tempDatabase);
    }

    /** Not instantiable. */
    private RamachandranDatabase()
    {
        throw new IllegalArgumentException("not instantiable");
    }

    private static class RamachandranUnit implements WorkUnit
    {
        public final String filename;
        public final Map<CustomKey,PreDistribution> targetMap;

        public RamachandranUnit(String filename, Map<CustomKey,PreDistribution> targetMap)
        {
            this.filename = filename;
            this.targetMap = targetMap;
        }

        public Result call()
        {
            // read data from zipped file
            try( GZIPInputStream gzip = new GZIPInputStream(new FileInputStream(filename));
                 BufferedReader br = new BufferedReader(new InputStreamReader(gzip));       )
                {
                    // keep track of the last line so we know if we've changed blocks
                    CustomAminoAcid lastCentralAminoAcid  = null;                // the amino acid in field 0
                    Direction lastDirection               = null;                // left or right in field 1
                    CustomAminoAcid lastAdjacentAminoAcid = null;                // the amino acid in field 2
                    
                    // temporary storage while reading a block
                    List<Double> tempPhis                = new LinkedList<>();  // backbone angle phi
                    List<Double> tempPsis                = new LinkedList<>();  // backbone angle psi
                    List<Double> tempProbabilities       = new LinkedList<>();  // log probabilities

                    while (true)
                        {
                            String currentLine = br.readLine();
                            
                            // break out when we have reached the end of the file
                            if ( currentLine == null )
                                break;
                            
                            // ignore comments and blank lines
                            String[] fields = currentLine.split("\\s+");
                            if ( currentLine.startsWith("#") || fields.length != 8)
                                continue;
                            
                            // special fix for proline -- PRO means trans proline
                            if ( fields[0].equals("PRO") )
                                fields[0] = "TPR";
                            if ( fields[2].equals("PRO") )
                                fields[2] = "TPR";

                            // parse to enum constants
                            CustomAminoAcid currentCentralAminoAcid = CustomAminoAcid.valueOf(fields[0]);
                            Direction currentDirection = Direction.valueOf(fields[1].toUpperCase());
                            CustomAminoAcid currentAdjacentAminoAcid = CustomAminoAcid.valueOf(fields[2]);

                            // detect a change in data block
                            if (  lastCentralAminoAcid != null &&
                                ( lastCentralAminoAcid  != currentCentralAminoAcid ||
                                  lastDirection         != currentDirection        ||
                                  lastAdjacentAminoAcid != currentAdjacentAminoAcid  ) )
                                {
                                    // create CustomKey object
                                    CustomKey customKey = new CustomKey(lastCentralAminoAcid, lastDirection, lastAdjacentAminoAcid);

                                    // create PreDistribution object
                                    short[] phiArray = Shorts.toArray(tempPhis);
                                    short[] psiArray = Shorts.toArray(tempPsis);
                                    float[] logProbabilityArray = Floats.toArray(tempProbabilities);
                                    PreDistribution preDistribution = new PreDistribution(phiArray, psiArray, logProbabilityArray);

                                    // add to database
                                    targetMap.put(customKey, preDistribution);

                                    // reset lists
                                    tempPhis = new LinkedList<>();
                                    tempPsis = new LinkedList<>();
                                    tempProbabilities = new LinkedList<>();
                                }

                            // parse fields and add to temporary lists
                            tempPhis.add(Double.valueOf(fields[3]));
                            tempPsis.add(Double.valueOf(fields[4]));
                            tempProbabilities.add(Double.valueOf(fields[6]));

                            // remember for the next line
                            lastCentralAminoAcid = currentCentralAminoAcid;
                            lastDirection = currentDirection;
                            lastAdjacentAminoAcid = currentAdjacentAminoAcid;
                        }

                    // deal with edge case
                    // create CustomKey object
                    CustomKey customKey = new CustomKey(lastCentralAminoAcid, lastDirection, lastAdjacentAminoAcid);

                    // create PreDistribution object
                    short[] phiArray = Shorts.toArray(tempPhis);
                    short[] psiArray = Shorts.toArray(tempPsis);
                    float[] logProbabilityArray = Floats.toArray(tempProbabilities);
                    PreDistribution preDistribution = new PreDistribution(phiArray, psiArray, logProbabilityArray);

                    // add to database
                    targetMap.put(customKey, preDistribution);
                }
            catch (IOException e)
                {
                    e.printStackTrace();
                    System.exit(1);
                }
            
            return null;
        }
    }

    /** Indicates whether the amino acid is to the left or to the right of the central one. */
    private enum Direction
    {
        /** amino acid is to the left of the central one */
        LEFT,

        /** amino acid is to the right of the central one */
        RIGHT;
    }

    /** Represents the amino acids in the Dunbrack library. */
    private enum CustomAminoAcid
    {
        ALL, ALA, GLY, VAL, LEU, ILE, TPR, CPR, PHE, TYR, TRP, SER, THR, CYS, MET, ASN, GLN, LYS, ARG, HIS, ASP, GLU;

        public static CustomAminoAcid getAminoAcid(AminoAcid aminoAcid, double omega)
        {
            if ( aminoAcid == AminoAcid.LPRO )
                {
                    if ( omega < 45.0 && omega > -45.0 )
                        return CPR;
                    else
                        return TPR;
                }
            else if ( aminoAcid == AminoAcid.DPRO )
                // treat D-proline as trans-L-proline
                return TPR;
            else if ( aminoAcid.rotamerType == AminoAcid.RotamerType.SPECIAL )
                throw new IllegalArgumentException("special amino acid");
            else
                return CustomAminoAcid.valueOf(aminoAcid.name());
        }

        public static CustomAminoAcid getAminoAcidNoCpr(AminoAcid aminoAcid, double omega)
        {
            CustomAminoAcid aa = getAminoAcid(aminoAcid, omega);
            if ( aa == CPR )
                return TPR;
            return aa;
        }
    }

    /**
     * Lightweight class for use as keys in the database hash table.  Everything is an enum so shouldn't
     * use much memory.
     */
    private static class CustomKey
    {
        /** the central amino acid whose (psi,phi) angles are the subject of the probability distribution */
        private CustomAminoAcid centralAminoAcid;

        /** whether the adjacent amino acid is to the left or right of the central one */
        private Direction direction;
        
        /** the adjacent amino acid whose identity will influence the (psi,phi) of the central residue */
        private CustomAminoAcid adjacentAminoAcid;

        /** simple constructor */
        public CustomKey(CustomAminoAcid centralAminoAcid, Direction direction, CustomAminoAcid adjacentAminoAcid)
        {
            this.centralAminoAcid = centralAminoAcid;
            this.direction = direction;
            this.adjacentAminoAcid = adjacentAminoAcid;
        }

        @Override
        public String toString()
        {
            return String.format("%5s %5s %5s\n", centralAminoAcid.toString(), direction.toString(), adjacentAminoAcid.toString() );
        }

        @Override
        public int hashCode()
        {
            return Objects.hash(centralAminoAcid, direction, adjacentAminoAcid);
        }

        @Override
        public boolean equals(Object obj)
        {
            if ( obj == null )
                return false;
            if ( obj == this )
                return true;
            if ( !(obj instanceof CustomKey) )
                return false;

            CustomKey another = (CustomKey)obj;
            if ( this.centralAminoAcid == another.centralAminoAcid &&
                 this.direction == another.direction &&
                 this.adjacentAminoAcid == another.adjacentAminoAcid )
                return true;
            return false;
        }
    }

    /**
     * Lightweight class that stores the neighbor-dependent Ramachandran data.  Can be converted to
     * DiscreteProbabilityDistribution.
     */
    private static class PreDistribution
    {
        /** the phi dihedral angles in degrees */
        private final short[] phis;

        /** the psi dihedral angles in degrees */
        private final short[] psis;

        /** log(probability) in no units */
        private final float[] logProbabilities;

        /** simple constructor */
        public PreDistribution(short[] phis, short[] psis, float[] logProbabilities)
        {
            this.phis = phis;
            this.psis = psis;
            this.logProbabilities = logProbabilities;
        }

        /**
         * Converts this PreDistribution to a DiscreteProbabilityDistribution:
         * outcomes: RotamerLibrary.Angles (phi,psi)
         * probabilities: are converted from log values to normal ones
         * @return the DiscreteProbabilityDistribution corresponding to this PreDistribution
         */
        public DiscreteProbabilityDistribution<RotamerLibrary.Angles> getDPD()
        {
            // convert primitive shorts to Doubles for use in the distribution
            List<RotamerLibrary.Angles> outcomes = new LinkedList<>();
            for (int i=0; i < phis.length; i++)
                {
                    Double thisPhi = Double.valueOf(phis[i]);
                    Double thisPsi = Double.valueOf(psis[i]);
                    RotamerLibrary.Angles theseAngles = new RotamerLibrary.Angles(thisPhi, thisPsi);
                    outcomes.add(theseAngles);
                }

            // turn log probabilities (floats) into normal probabilities (Doubles)
            List<Double> probabilities = new LinkedList<>();
            for (float f : logProbabilities)
                {
                    Double thisLogValue = -1.0 * Double.valueOf(f);
                    Double thisProbability = Math.exp(thisLogValue);
                    probabilities.add(thisProbability);
                }

            return new DiscreteProbabilityDistribution<RotamerLibrary.Angles>(outcomes, probabilities);
        }

        @Override
        public String toString()
        {
            String returnString = "[";
            int n = phis.length;
            for (int i=0; i < n - 1; i++)
                returnString = returnString + String.format("%5d %5d %.6f,\n", phis[i], psis[i], logProbabilities[i]);
            returnString = returnString + String.format("%5d %5d %.6f]", phis[n-1], psis[n-1], logProbabilities[n-1]);
            return returnString;
        }
        
        @Override
        public int hashCode()
        {
            return Objects.hash(phis,psis,logProbabilities);
        }

        @Override
        public boolean equals(Object obj)
        {
            if ( obj == null )
                return false;
            if ( obj == this )
                return true;
            if ( !(obj instanceof PreDistribution) )
                return false;

            PreDistribution another = (PreDistribution)obj;
            if ( Arrays.equals(phis, another.phis) &&
                 Arrays.equals(psis, another.psis) &&
                 Arrays.equals(logProbabilities, another.logProbabilities) )
                return true;
            return false;
        }
    }

    /**
     * Finds an entry in the database.  Result is a PreDistribution, which holds the backbone angles (phi,psi) as
     * shorts and the log of probability as a float.  This can be converted to a DiscreteProbabilityDistribution
     * on demand.
     * @param centralAminoAcid the amino acid in the middle
     * @param centralOmega the omega dihedral angle of the amino acid in the middle in degrees
     * @param direction whether the adjacent amino acid is to the left or right of the central one
     * @param adjacentAminoAcid the amino acid next to the central one
     * @param adjacentOmega the omega dihedral angle of the adjacent amino acid
     * @return the conditional probability for (phi,psi) given this pair of amino acids
     */
    private static PreDistribution locate(AminoAcid centralAminoAcid, double centralOmega, Direction direction, AminoAcid adjacentAminoAcid, double adjacentOmega)
    {
        CustomAminoAcid centralAA = CustomAminoAcid.getAminoAcid(centralAminoAcid, centralOmega);
        CustomAminoAcid adjacentAA = CustomAminoAcid.getAminoAcidNoCpr(adjacentAminoAcid, adjacentOmega);
        CustomKey thisKey = new CustomKey(centralAA, direction, adjacentAA);
        PreDistribution result = DATABASE.get(thisKey);
        if ( result == null )
            {
                String errorString = String.format("central: %s  direction: %s  adjacent: %s ", centralAminoAcid.toString(), direction.toString(), adjacentAminoAcid.toString());
                errorString += String.format(" central: %s  adjacent: %s ", centralAA.toString(), adjacentAA.toString());
                throw new NullPointerException("Could not locate the requsted target (" + errorString + ")!");
            }
        return result;
    }

    /**
     * Finds an entry in the database.  Result is a PreDistribution, which holds the backbone angles (phi,psi) as
     * shorts and the log of probability as a float.  This can be converted to a DiscreteProbabilityDistribution
     * on demand.
     *
     * This version finds Xxx-ALL or ALL-Xxx conditional probabilities.
     * @param centralAminoAcid the amino acid in the middle
     * @param centralOmega the omega dihedral angle of the amino acid in the middle in degrees
     * @param direction whether the adjacent amino acid is to the left or right of the central one
     * @return the conditional probability for (phi,psi) given this amino acid and anything to the left or right
     */
    private static PreDistribution locateAll(AminoAcid centralAminoAcid, double centralOmega, Direction direction)
    {
        CustomAminoAcid centralAA = CustomAminoAcid.getAminoAcid(centralAminoAcid, centralOmega);
        CustomKey thisKey = new CustomKey(centralAA, direction, CustomAminoAcid.ALL);
        PreDistribution result = DATABASE.get(thisKey);
        if ( result == null )
            {
                String errorString = String.format("central: %s  direction: %s  adjacent: ALL", centralAminoAcid.toString(), direction.toString());
                throw new NullPointerException("Could not locate the requsted target (" + errorString + ")!");
            }
        return result;
    }

    /**
     * Returns a DiscreteProbabilityDistribution Pr(phi,psi | left AA - central AA).  That is, given a sequence (left AA, central AA)
     * in the N to C direction, gives the conditional probability that the central AA will have backbone angles (phi,psi). <p>
     * Parameters are not checked for correctness!  For example, asking for the central amino acid to be AminoAcid.ALL will result
     * in a NullPointerException.  Probability distributions are generated on demand from internal lists of primitives to save memory.
     * @param leftAminoAcid the amino acid on the left in the N to C direction
     * @param leftOmega the omega dihedral angle in degrees for the left amino acid
     * @param centralAminoAcid the amino acid on the right in the N to C direction
     * @param centralOmega the omega dihedral angle in degrees for the central amino acid
     * @return the DiscreteProbabilityDistribution of (psi,phi) values of the central amino acid
     */
    public static DiscreteProbabilityDistribution<RotamerLibrary.Angles> getLeftDistribution(AminoAcid leftAminoAcid, double leftOmega,
                                                                                             AminoAcid centralAminoAcid, double centralOmega)
    {
        PreDistribution preDistribution = locate(centralAminoAcid, centralOmega, Direction.LEFT, leftAminoAcid, leftOmega);
        return preDistribution.getDPD();
    }

    /**
     * Returns a DiscreteProbabilityDistribution Pr(phi,psi | central AA - right AA).  That is, given a sequence (central AA, right AA)
     * in the N to C direction, gives the conditional probability that the central AA will have backbone angles (phi,psi). <p>
     * Parameters are not checked for correctness!  For example, asking for the central amino acid to be AminoAcid.ALL will result
     * in a NullPointerException.  Probability distributions are generated on demand from internal lists of primitives to save memory.
     * @param centralAminoAcid the amino acid on the left in the N to C direction
     * @param centralOmega the omega dihedral angle in degrees for the central amino acid
     * @param rightAminoAcid the amino acid on the right in the N to C direction
     * @param rightOmega the omega dihedral angle in degrees for the right amino acid
     * @return the DiscreteProbabilityDistribution of (psi,phi) values of the central amino acid
     */
    public static DiscreteProbabilityDistribution<RotamerLibrary.Angles> getRightDistribution(AminoAcid centralAminoAcid, double centralOmega,
                                                                                              AminoAcid rightAminoAcid, double rightOmega)
    {
        PreDistribution preDistribution = locate(centralAminoAcid, centralOmega, Direction.RIGHT, rightAminoAcid, rightOmega);
        return preDistribution.getDPD();
    }

    /**
     * Returns a DiscreteProbabilityDistribution Pr(phi,psi | left AA - central AA - right AA).  That is, given
     * a sequence (left AA, central AA, right AA) in the N to C direction, gives the conditional probability that
     * the central AA will have backbone angles (phi,psi).<p>
     * <p>
     * This assumes all the amino acids are in the database.  Don't use AminoAcid.ALL!  This is not checked and
     * will throw a NullPointerException.<p>
     * <p>
     * Algorithm, which is executed on demand:<p>
     *  1. For all (psi,phi), calculate log Pr(phi,psi|C,L) + log Pr(phi,psi|C,R) - log Pr(phi,psi|C,R=ALL).<p>
     *  2. Convert back to a regular probability.<p>
     *  3. Create a DiscreteProbabilityDistribution, which will normalize the probabilities automatically.<p>
     * <p>
     * Pr(phi,psi|C,L) means the conditional probability that (phi,psi) will be observed given that the central
     * amino acid is C and an amino acid L is to the left of it.  The algorithm amounts to multiplying the independent
     * probabilities that we get C,L and C,R, conditional on the central residue being C.
     *
     * @param leftAminoAcid the amino acid on the left on the left in the N to C direction
     * @param leftOmega the omega dihedral angle in degrees for the left amino acid
     * @param centralAminoAcid the amino acid in the middle in the N to C direction
     * @param centralOmega the omega dihedral angle in degrees for the central amino acid
     * @param rightAminoAcid the amino acid on the right in the N to C direction
     * @param rightOmega the omega dihedral angle in degrees for the right amino acid
     * @return the DiscreteProbabilityDistribution of (psi,phi) values of the right amino acid
     */
    public static DiscreteProbabilityDistribution<RotamerLibrary.Angles>
                    getTripletDistribution(AminoAcid leftAminoAcid, double leftOmega,
                                           AminoAcid centralAminoAcid, double centralOmega,
                                           AminoAcid rightAminoAcid, double rightOmega)
    {
        // get the relevant data
        PreDistribution leftDistribution  =    locate(centralAminoAcid, centralOmega, Direction.LEFT,  leftAminoAcid,  leftOmega);
        PreDistribution rightDistribution =    locate(centralAminoAcid, centralOmega, Direction.RIGHT, rightAminoAcid, centralOmega);
        PreDistribution allDistribution   = locateAll(centralAminoAcid, centralOmega, Direction.RIGHT);

        // calculate new log sums and convert to regular probabilities
        int n = leftDistribution.phis.length;
        List<Double> newProbabilities = new LinkedList<>();
        double sum = 0.0;
        for (int i=0; i < n; i++)
            {
                double logProbabilitySum = (double)leftDistribution.logProbabilities[i] +
                                           (double)rightDistribution.logProbabilities[i] -
                                           (double)allDistribution.logProbabilities[i];
                Double newProbability = Math.exp(-1.0*logProbabilitySum);
                sum += newProbability;
                newProbabilities.add(newProbability);
            }

        // we're going to normalize anyways to make the debugging easier
        List<Double> normalizedProbabilities = new LinkedList<>();
        for (Double d : newProbabilities)
            normalizedProbabilities.add(d/sum);

        // convert primitive shorts to Doubles for use in the distribution
        List<RotamerLibrary.Angles> outcomes = new LinkedList<>();
        for (int i=0; i < n; i++)
            {
                Double thisPhi = Double.valueOf(leftDistribution.phis[i]);
                Double thisPsi = Double.valueOf(leftDistribution.psis[i]);
                RotamerLibrary.Angles theseAngles = new RotamerLibrary.Angles(thisPhi, thisPsi);
                outcomes.add(theseAngles);
            }

        // create DiscreteProbabilityDistribution
        return new DiscreteProbabilityDistribution<RotamerLibrary.Angles>(outcomes, normalizedProbabilities);
    }

    /** Forces the static initializer to run. */
    public static void load()
    {
        System.out.printf("RamachandranDatabase loaded with %d entries.\n", DATABASE.size());
    }

    /** For testing. */
    public static void main(String[] args)
    {
        System.out.println("\n" + DATABASE.size());

        // iterate through all possibilities 
        List<AminoAcid> aminoAcids = new ArrayList<>();
        for (AminoAcid a : AminoAcid.values())
            aminoAcids.add(a);
        aminoAcids.remove(AminoAcid.DPRO);
        aminoAcids.remove(AminoAcid.TS);
        for (int i=0; i < aminoAcids.size(); i++)
            {
                AminoAcid aa1 = aminoAcids.get(i);
                for (int j=0; j < aminoAcids.size(); j++)
                    {
                        AminoAcid aa2 = aminoAcids.get(j);
                        List<Double> list1 = ImmutableList.of(180.0);
                        List<Double> list2 = ImmutableList.of(180.0);
                        if ( aa1 == AminoAcid.LPRO )
                            list1 = ImmutableList.of(180.0, 0.0);
                        if ( aa2 == AminoAcid.LPRO )
                            list2 = ImmutableList.of(180.0, 0.0);
                        for (Double dihedral1 : list1)
                            {
                                for (Double dihedral2 : list2)
                                    {
                                        System.out.printf("%5s (%3.0f) - LEFT  - %5s (%3.0f)\n", aa1.shortName, dihedral1, aa2.shortName, dihedral2);
                                        PreDistribution d = locate(aa1, dihedral1, Direction.LEFT, aa2, dihedral2);
                                        System.out.printf("%5s (%3.0f) - RIGHT - %5s (%3.0f)\n", aa1.shortName, dihedral1, aa2.shortName, dihedral2);
                                        PreDistribution d2 = locate(aa1, dihedral1, Direction.RIGHT, aa2, dihedral2);
                                    }
                            }
                    }
            }

        /*//System.out.println(locate(AminoAcid.ALA, Direction.LEFT, AminoAcid.ALL));
        //System.out.println(locate(AminoAcid.ALA, Direction.LEFT, AminoAcid.ALL).getDPD());
        //System.out.println(getLeftDistribution(AminoAcid.ALA, AminoAcid.ALL));
        //System.out.println(getRightDistribution(AminoAcid.ALA, AminoAcid.ALL));

	    //System.out.println(locate(AminoAcid.VAL,Direction.RIGHT,AminoAcid.VAL));

        // note: probabilities will not appear to add to 1.0 even though they actually do because
        // not all entries in the DPD will be printed out!
        
        for (int i=0; i<1000 ; i++)
        {
            DiscreteProbabilityDistribution<RotamerLibrary.Angles> dpd = getTripletDistribution(AminoAcid.LEU, AminoAcid.ASN, AminoAcid.CPR);
            System.out.println("Size of dpd before: " + dpd.getSize());
            RotamerLibrary.Angles newPhiPsi = getNearbyPhiPsi(new RotamerLibrary.Angles(50.0,50.0), 100.0, dpd);
            System.out.println(newPhiPsi);
        }
        */
/*
        // phi,psi for gly in pro-gly
        DiscreteProbabilityDistribution<RotamerLibrary.Angles> dpdGly = getLeftDistribution(AminoAcid.PRO,AminoAcid.GLY);
        String description = "";
        double sum = 0.0;
        for (int i=0; i < dpdGly.outcomes.size(); i++)
            {
                Double probability = dpdGly.inputProbabilities.get(i) * 100.0;
                RotamerLibrary.Angles angles = dpdGly.outcomes.get(i);
                double phi = -1.0 * angles.getPhi();
                double psi = -1.0 * angles.getPsi();
                description += String.format("%6.1f,%6.1f,%9.6f\n", phi, psi, probability);
            }
        InputFileFormat.writeStringToDisk(description,"gly.csv");

        // phi,psi for pro in pro-gly
        DiscreteProbabilityDistribution<RotamerLibrary.Angles> dpdPro = getRightDistribution(AminoAcid.PRO,AminoAcid.GLY);
        description = "";
        for (int i=0; i < dpdPro.outcomes.size(); i++)
            {
                Double probability = dpdPro.inputProbabilities.get(i) * 100.0;
                RotamerLibrary.Angles angles = dpdPro.outcomes.get(i);
                double phi = -1.0 * angles.getPhi();
                double psi = -1.0 * angles.getPsi();
                description += String.format("%6.1f,%6.1f,%9.6f\n", phi, psi, probability);
                if ( phi >= 51.0 && phi <= 108.0 && psi >= -10.0 && psi <= 47.0 )
                sum += probability;
            }
        System.out.println(sum);
        InputFileFormat.writeStringToDisk(description,"pro.csv");

        

        
        //System.out.println(getTripletDistribution(AminoAcid.TYR, AminoAcid.ARG, AminoAcid.TRP).toDebugString(0.0001,10));
        //System.out.println(getTripletDistribution(AminoAcid.PHE, AminoAcid.ARG, AminoAcid.TRP).toDebugString(0.0001,10));
        
        //System.out.println(getTripletDistribution(AminoAcid.ASN, AminoAcid.TPR, AminoAcid.SER).toDebugString(0.0001,10));
        //System.out.println(getTripletDistribution(AminoAcid.ASN, AminoAcid.TPR, AminoAcid.THR).toDebugString(0.0001,10));
        Scanner scanner = new Scanner(System.in);
        System.out.println("Press enter to continue.");
        scanner.nextLine();
  */      
    }
}
