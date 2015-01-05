import java.util.*;
import java.util.zip.GZIPInputStream;
import java.io.*;
import com.google.common.collect.*;
import com.google.common.primitives.*;

/**
 * This singleton stores data on how (phi,psi) depends on the neighboring residue identity based on the Dunbrack data.<p>
 * Reference: Daniel Ting, Guoli Wang, Maxim Shapovalov, Rajib Mitra, Michael I. Jordan, Roland L. Dunbrack, Jr.
 * Neighbor-dependent Ramachandran probability distributions of amino acids developed from a<u>
 * hierarchical Dirichlet process model.</u><i> PLOS Comp. Biol.</i><em> (April 2010)</em>.<p><p>
 * Use this command to check memory usage:<p>
 * <code>jps | grep Rama | awk '{print $1}' | xargs jmap -histo:live | awk '{if ($1 == "Total" || NR < 20) {print}}'
 */
public class RamachandranLibrary implements Singleton
{
    /** the singleton instance */
    public static final RamachandranLibrary INSTANCE = new RamachandranLibrary();
    
    /**
     * stores the Ramachandran data
     */
    private final Map<CustomKey,PreDistribution> database;

    /** read all the Ramachandran data */
    private RamachandranLibrary()
    {
        if (INSTANCE != null)
            throw new IllegalStateException("this should be a singleton!");
        
        // temporary copy of the database
        Map<CustomKey,PreDistribution> tempDatabase = new HashMap<>();
        
        // read data from zipped file
        GZIPInputStream gzip = null;
        BufferedReader br = null;
        try
            {
                gzip = new GZIPInputStream(new FileInputStream(Settings.RAMACHANDRAN_DATA_FILENAME));
                br = new BufferedReader(new InputStreamReader(gzip));

                // keep track of the last line so we know if we've changed blocks
                LibraryValue lastCentralAminoAcid  = LibraryValue.ALA;        // the amino acid in field 0
                Direction lastDirection         = Direction.LEFT;       // left or right in field 1
                LibraryValue lastAdjacentAminoAcid = LibraryValue.ALL;        // the amino acid in field 2
                
                // temporary storage while reading a block
                List<Double> tempPhis                = new LinkedList<>();  // backbone angle phi
                List<Double> tempPsis                = new LinkedList<>();  // backbone angle psi
                List<Double> tempProbabilities       = new LinkedList<>();  // log probabilities

                System.out.println("Reading Ramachandran database...");
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
                        LibraryValue currentCentralAminoAcid = LibraryValue.valueOf(fields[0]);
                        Direction currentDirection = Direction.valueOf(fields[1].toUpperCase());
                        LibraryValue currentAdjacentAminoAcid = LibraryValue.valueOf(fields[2]);

                        // for debugging only -- shortens runtime
                        //if ( ! fields[0].equals("ALA") )
                        //    break;

                        // detect a change in data block
                        if ( lastCentralAminoAcid  != currentCentralAminoAcid ||
                             lastDirection         != currentDirection        ||
                             lastAdjacentAminoAcid != currentAdjacentAminoAcid )
                            {
                                // print status
                                //System.out.print( String.format("%5s %5s %5s\r", lastCentralAminoAcid.toString(),
                                //                                                 lastDirection.toString(),
                                //                                                 lastAdjacentAminoAcid.toString() ) );

                                // create CustomKey object
                                CustomKey customKey = new CustomKey(lastCentralAminoAcid, lastDirection, lastAdjacentAminoAcid);

                                // create PreDistribution object
                                short[] phiArray = Shorts.toArray(tempPhis);
                                short[] psiArray = Shorts.toArray(tempPsis);
                                float[] logProbabilityArray = Floats.toArray(tempProbabilities);
                                PreDistribution preDistribution = new PreDistribution(phiArray, psiArray, logProbabilityArray);

                                // add to database
                                tempDatabase.put(customKey, preDistribution);

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
                tempDatabase.put(customKey, preDistribution);

                br.close();
                gzip.close();
            }
        catch (IOException e)
            {
                e.printStackTrace();
                System.exit(1);
            }

        // make dataset immutable
        /*for (CustomKey k : tempDatabase.keySet())
        {
            String description = k.toString() + " / " + k.hashCode();;
            if ( description.indexOf("Pro") > -1 && description.indexOf("Asn") > -1 )
                System.out.println(description);
        }*/
        database = ImmutableMap.copyOf(tempDatabase);
    }

    /** Indicates whether the amino acid is to the left or to the right of the central one. */
    private enum Direction
    {
        /** amino acid is to the left of the central one */
        LEFT,
        /** amino acid is to the right of the central one */
        RIGHT;
    }

    /** Indicates the possible library values for an amino acid. Can include ALL to signify all amino acids. */
    private enum LibraryValue
    {
        ALA("Ala"),
        GLY("Gly"),
        VAL("Val"),
        LEU("Leu"),
        ILE("Ile"),

        PRO("Pro"),
        CPR("Cpr"),
        TPR("Tpr"),
        DPRO("DPro"),

        PHE("Phe"),
        TYR("Tyr"),
        TRP("Trp"),
        SER("Ser"),
        THR("Thr"),
        CYS("Cys"),
        MET("Met"),
        ASN("Asn"),
        GLN("Gln"),
        LYS("Lys"),
        ARG("Arg"),
        HIS("His"),
        ASP("Asp"),
        GLU("Glu"),
        ALL("All");

        /** String containing the amino acid name in the AminoAcid enum */
        public final String aminoAcidName;

        LibraryValue(String shortName)
        {
            this.aminoAcidName = shortName;
        }

        public AminoAcid getAminoAcid()
        {
            if (this.equals(LibraryValue.PRO) || this.equals(LibraryValue.ALL))
                throw new IllegalArgumentException("Not part of amino acid list");

            return AminoAcid.valueOf(aminoAcidName.toUpperCase());
        }

        public static LibraryValue getLibraryValue(AminoAcid aminoAcid)
        {
            return LibraryValue.valueOf(aminoAcid.getCompatible().name());
        }

        public String toString()
        {
            return aminoAcidName;
        }
    
    }

    /**
     * Lightweight class for use as keys in the database hash table.  Everything is an enum so shouldn't
     * use much memory.
     */
    private static class CustomKey
    {
        /** the central amino acid whose (psi,phi) angles are the subject of the probability distribution */
        private LibraryValue centralAminoAcid;

        /** whether the adjacent amino acid is to the left or right of the central one */
        private Direction direction;
        
        /** the adjacent amino acid whose identity will influence the (psi,phi) of the central residue */
        private LibraryValue adjacentAminoAcid;

        /** simple constructor */
        public CustomKey(LibraryValue centralAminoAcid, Direction direction, LibraryValue adjacentAminoAcid)
        {
            this.centralAminoAcid = centralAminoAcid;
            this.direction = direction;
            this.adjacentAminoAcid = adjacentAminoAcid;
        }

        /**
         * Give a simple text representation of this key.
         * @return the description
         */
        public String toString()
        {
            return String.format("%5s %5s %5s\n", centralAminoAcid.toString(),
                                                       direction.toString(),
                                                       adjacentAminoAcid.toString() );
        }

        /**
         * Returns the hash code.
         * @return the hash code
         */
        public int hashCode()
        {
            return Objects.hash(centralAminoAcid, direction, adjacentAminoAcid);
        }

        /**
         * Simple test for logical equivalence.
         * @param obj the object we are comparing to
         * @return whether the objects are equal
         */
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
         * outcomes: SidechainRotamerLibrary.BackboneAngles (phi,psi)
         * probabilities: are converted from log values to normal ones
         * @return the DiscreteProbabilityDistribution corresponding to this PreDistribution
         */
        public DiscreteProbabilityDistribution<SidechainRotamerLibrary.BackboneAngles> getDPD()
        {
            // convert primitive shorts to Doubles for use in the distribution
            List<SidechainRotamerLibrary.BackboneAngles> outcomes = new LinkedList<>();
            for (int i=0; i < phis.length; i++)
                {
                    Double thisPhi = Double.valueOf(phis[i]);
                    Double thisPsi = Double.valueOf(psis[i]);
                    SidechainRotamerLibrary.BackboneAngles theseAngles = new SidechainRotamerLibrary.BackboneAngles(thisPhi, thisPsi);
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

            return new DiscreteProbabilityDistribution<SidechainRotamerLibrary.BackboneAngles>(outcomes, probabilities);
        }

        /**
         * Give a simple text representation of this key.
         * @return the description
         */
        public String toString()
        {
            String returnString = "[";
            int n = phis.length;
            for (int i=0; i < n - 1; i++)
                returnString = returnString + String.format("%5d %5d %.6f,\n", phis[i], psis[i], logProbabilities[i]);
            returnString = returnString + String.format("%5d %5d %.6f]", phis[n-1], psis[n-1], logProbabilities[n-1]);
            return returnString;
        }

        /**
         * Returns the hash code.
         * @return the hash code
         */
        public int hashCode()
        {
            return Objects.hash(phis,psis,logProbabilities);
        }

        /**
         * Simple test for logical equivalence.
         * @param obj the object we are comparing to
         * @return whether the objects are equal
         */
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
     * Gives a short description of this singleton.
     * @return the short description
     */
    public String toString()
    {
        return String.format("Ramachandran database loaded with %d entries", database.size());
    }

    /**
     * Finds an entry in the database.  Result is a PreDistribution, which holds the backbone angles (phi,psi) as
     * shorts and the log of probability as a float.  This can be converted to a DiscreteProbabilityDistribution
     * on demand.
     * @param centralAminoAcid the amino acid in the middle
     * @param direction whether the adjacent amino acid is to the left or right of the central one
     * @param adjacentAminoAcid the amino acid next to the central one
     * @return the conditional probability for (phi,psi) given this pair of amino acids
     */
    private static PreDistribution locate(LibraryValue centralAminoAcid, Direction direction, LibraryValue adjacentAminoAcid)
    {
        CustomKey thisKey = new CustomKey(centralAminoAcid, direction, adjacentAminoAcid);
        PreDistribution result = INSTANCE.database.get(thisKey);
        if ( result == null )
            {
                // there is no information for adjacent cis-proline
                if ( adjacentAminoAcid == LibraryValue.CPR )
                    {
                        thisKey = new CustomKey(centralAminoAcid, direction, AminoAcid.TPR);
                        result = INSTANCE.database.get(thisKey);
                    }

                if ( result == null )
                    {
                        String errorString = String.format("central: %s  direction: %s  adjacent: %s", centralAminoAcid.toString(), direction.toString(), adjacentAminoAcid.toString());
                        throw new NullPointerException("Could not locate the requsted target (" + errorString + ")!");
                    }
            }
        return result;
    }

    /**
     * Returns a DiscreteProbabilityDistribution Pr(phi,psi | left AA - central AA).  That is, given a sequence (left AA, central AA)
     * in the N to C direction, gives the conditional probability that the central AA will have backbone angles (phi,psi). <p>
     * Parameters are not checked for correctness!  For example, asking for the central amino acid to be AminoAcid.ALL will result
     * in a NullPointerException.  Probability distributions are generated on demand from internal lists of primitives to save memory.
     * @param leftAminoAcid the amino acid on the left in the N to C direction
     * @param centralAminoAcid the amino acid on the right in the N to C direction
     * @return the DiscreteProbabilityDistribution of (psi,phi) values of the central amino acid
     */
    public static DiscreteProbabilityDistribution<SidechainRotamerLibrary.BackboneAngles> getLeftDistribution(AminoAcid leftAminoAcid, AminoAcid centralAminoAcid)
    {
        PreDistribution preDistribution = locate(LibraryValue.getLibraryValue(centralAminoAcid), Direction.LEFT, LibraryValue.getLibraryValue(leftAminoAcid));
        return preDistribution.getDPD();
    }

    /**
     * Returns a DiscreteProbabilityDistribution Pr(phi,psi | central AA - right AA).  That is, given a sequence (central AA, right AA)
     * in the N to C direction, gives the conditional probability that the central AA will have backbone angles (phi,psi). <p>
     * Parameters are not checked for correctness!  For example, asking for the central amino acid to be AminoAcid.ALL will result
     * in a NullPointerException.  Probability distributions are generated on demand from internal lists of primitives to save memory.
     * @param centralAminoAcid the amino acid on the left in the N to C direction
     * @param rightAminoAcid the amino acid on the right in the N to C direction
     * @return the DiscreteProbabilityDistribution of (psi,phi) values of the central amino acid
     */
    public static DiscreteProbabilityDistribution<SidechainRotamerLibrary.BackboneAngles> getRightDistribution(AminoAcid centralAminoAcid, AminoAcid rightAminoAcid)
    {
        PreDistribution preDistribution = locate(LibraryValue.getLibraryValue(centralAminoAcid), Direction.RIGHT, LibraryValue.getLibraryValue(rightAminoAcid));
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
     * @param centralAminoAcid the amino acid in the middle in the N to C direction
     * @param rightAminoAcid the amino acid on the right in the N to C direction
     * @return the DiscreteProbabilityDistribution of (psi,phi) values of the central amino acid
     */
    public static DiscreteProbabilityDistribution<SidechainRotamerLibrary.BackboneAngles>
                    getTripletDistribution(AminoAcid leftAminoAcid, AminoAcid centralAminoAcid, AminoAcid rightAminoAcid)
    {
        // get the relevant data
        PreDistribution leftDistribution  = locate(LibraryValue.getLibraryValue(centralAminoAcid), Direction.LEFT, LibraryValue.getLibraryValue(leftAminoAcid));
        PreDistribution rightDistribution = locate(LibraryValue.getLibraryValue(centralAminoAcid), Direction.RIGHT, LibraryValue.getLibraryValue(rightAminoAcid));
        PreDistribution allDistribution   = locate(LibraryValue.getLibraryValue(centralAminoAcid), Direction.RIGHT, LibraryValue.ALL);

        //System.out.println(leftDistribution.phis.length + " " + rightDistribution.phis.length + " " + allDistribution.phis.length);
        //System.out.println(leftDistribution.psis.length + " " + rightDistribution.psis.length + " " + allDistribution.psis.length);
        //System.out.println(leftDistribution.logProbabilities.length + " " + rightDistribution.logProbabilities.length + " " + allDistribution.logProbabilities.length);

        //for (int i=0; i < leftDistribution.phis.length; i++)
        //    {
        //        System.out.println(leftDistribution.phis[i] + " " + rightDistribution.phis[i] + " " + allDistribution.phis[i]);
        //        System.out.println(leftDistribution.psis[i] + " " + rightDistribution.psis[i] + " " + allDistribution.psis[i]);
        //    }

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
        List<SidechainRotamerLibrary.BackboneAngles> outcomes = new LinkedList<>();
        for (int i=0; i < n; i++)
            {
                Double thisPhi = Double.valueOf(leftDistribution.phis[i]);
                Double thisPsi = Double.valueOf(leftDistribution.psis[i]);
                SidechainRotamerLibrary.BackboneAngles theseAngles = new SidechainRotamerLibrary.BackboneAngles(thisPhi, thisPsi);
                outcomes.add(theseAngles);
            }

        /*int count = 0;
        for (Double d : normalizedProbabilities)
            {
                count++;
                System.out.print(String.format("%.6f  ",d));
                if ( count == 20 )
                    {
                        count = 0;
                        System.out.println();
                    }
            }
        System.exit(1);*/

        // create DiscreteProbabilityDistribution
        return new DiscreteProbabilityDistribution<SidechainRotamerLibrary.BackboneAngles>(outcomes, normalizedProbabilities);
    }

    /** A method that returns nearby phi, psi by refactoring the Discrete Probability Distribution based on an input distance
    *   @param currentPhiPsi the SidechainRotamerLibrary.BackboneAngles that is used to find close phi-psi combinations for tweaking
    *   @param distance the Euclidean distance capturing the difference between two phi-psi pairs that is allowed for tweaking
    *   @param dpd the Discrete Probability Distribution that contains all the values for the phi-psi in the Ramachandran Library
    *   @return a new set of phi-psi angles */
    public static SidechainRotamerLibrary.BackboneAngles getNearbyPhiPsi(SidechainRotamerLibrary.BackboneAngles currentPhiPsi, double distance, DiscreteProbabilityDistribution<SidechainRotamerLibrary.BackboneAngles> dpd)
    {
        //initialize variables neeeded to create new discrete probability distribution
        List<SidechainRotamerLibrary.BackboneAngles> newOutcomes = new LinkedList<>();
        List<Double> probabilities = new LinkedList<>();
        
        //Refactor discrete probability distribution based on distance
        for (SidechainRotamerLibrary.BackboneAngles phiPsi : dpd.outcomes)
        {
            double diff = Math.sqrt(Math.pow(currentPhiPsi.getPhi()-phiPsi.getPhi(),2)+Math.pow(currentPhiPsi.getPsi()-phiPsi.getPsi(),2));
            if (diff < distance)
            {
                newOutcomes.add(phiPsi);
                double currentProbability = dpd.inputProbabilities.get(dpd.outcomes.indexOf(phiPsi));
                probabilities.add(currentProbability);
            }
        }

        //Create a new discrete probability distribution
        DiscreteProbabilityDistribution<SidechainRotamerLibrary.BackboneAngles> tweakDPD = new DiscreteProbabilityDistribution<>(newOutcomes, probabilities);

        //System.out.println(tweakDPD);
        //System.out.println("New DPD size: " + tweakDPD.getSize());

        //Draw randomly from discrete probability distribution
        SidechainRotamerLibrary.BackboneAngles newPhiPsi = tweakDPD.getRandom();
        return newPhiPsi;
    }

    /**
     * Returns the (phi,psi) distribution for a residue based on the identity of its neighbors.
     * @param peptide the peptide to analyze
     * @param residueIndex the index of the residue to get the (phi,psi) distribution of
     * @return the (phi,psi) distribution
     */
    public static DiscreteProbabilityDistribution<SidechainRotamerLibrary.BackboneAngles> getDistribution(Peptide peptide, int residueIndex)
    {
        // get adjacent amino acids
        Residue residue = peptide.sequence.get(residueIndex);
        List<Residue> sequence = peptide.sequence;
        AminoAcid centralAminoAcid = residue.aminoAcid;

        // no data are available for D-proline
        if ( centralAminoAcid == AminoAcid.DPRO )
            throw new IllegalArgumentException("doesn't handle D-proline");

        // for normal amino acids
        AminoAcid leftAminoAcid = null;
        if ( residueIndex > 0 )
            leftAminoAcid = sequence.get(residueIndex-1).aminoAcid;

        AminoAcid rightAminoAcid = null;
        if ( residueIndex < sequence.size() - 1 )
            rightAminoAcid = sequence.get(residueIndex+1).aminoAcid;

        // if D-proline is adjacent, set to null because we don't have data for that
        if ( leftAminoAcid == AminoAcid.DPRO )
            leftAminoAcid = null;
        if ( rightAminoAcid == AminoAcid.DPRO )
            rightAminoAcid = null;

        // obtain appropriate probability distribution of phi,psi
        DiscreteProbabilityDistribution<SidechainRotamerLibrary.BackboneAngles> DPD = null;
        if ( leftAminoAcid != null && rightAminoAcid != null )
            DPD = getTripletDistribution(leftAminoAcid, centralAminoAcid, rightAminoAcid);
        else if ( leftAminoAcid == null && rightAminoAcid != null )
            DPD = getRightDistribution(centralAminoAcid, rightAminoAcid);
        else if ( leftAminoAcid != null && rightAminoAcid == null )
            DPD = getLeftDistribution(leftAminoAcid, centralAminoAcid);
        else
            throw new IllegalArgumentException("should not be possible");
    
        // return result
        return DPD;
    }

    /** For testing. */
    public static void main(String[] args)
    {
        System.out.println("\nOperation complete.");
        //System.out.println(locate(AminoAcid.ALA, Direction.LEFT, AminoAcid.ALL));
        //System.out.println(locate(AminoAcid.ALA, Direction.LEFT, AminoAcid.ALL).getDPD());
        //System.out.println(getLeftDistribution(AminoAcid.ALA, AminoAcid.ALL));
        //System.out.println(getRightDistribution(AminoAcid.ALA, AminoAcid.ALL));

	    //System.out.println(locate(AminoAcid.VAL,Direction.RIGHT,AminoAcid.VAL));

        // note: probabilities will not appear to add to 1.0 even though they actually do because
        // not all entries in the DPD will be printed out!
        /*
        for (int i=0; i<1000 ; i++)
        {
            DiscreteProbabilityDistribution<SidechainRotamerLibrary.BackboneAngles> dpd = getTripletDistribution(AminoAcid.LEU, AminoAcid.ASN, AminoAcid.CPR);
            System.out.println("Size of dpd before: " + dpd.getSize());
            SidechainRotamerLibrary.BackboneAngles newPhiPsi = getNearbyPhiPsi(new SidechainRotamerLibrary.BackboneAngles(50.0,50.0), 100.0, dpd);
            System.out.println(newPhiPsi);
        }
        */

        // phi,psi for gly in pro-gly
        DiscreteProbabilityDistribution<SidechainRotamerLibrary.BackboneAngles> dpdGly = getLeftDistribution(AminoAcid.PRO,AminoAcid.GLY);
        String description = "";
        double sum = 0.0;
        for (int i=0; i < dpdGly.outcomes.size(); i++)
            {
                Double probability = dpdGly.inputProbabilities.get(i) * 100.0;
                SidechainRotamerLibrary.BackboneAngles angles = dpdGly.outcomes.get(i);
                double phi = -1.0 * angles.getPhi();
                double psi = -1.0 * angles.getPsi();
                description += String.format("%6.1f,%6.1f,%9.6f\n", phi, psi, probability);
            }
        InputFileFormat.writeStringToDisk(description,"gly.csv");

        // phi,psi for pro in pro-gly
        DiscreteProbabilityDistribution<SidechainRotamerLibrary.BackboneAngles> dpdPro = getRightDistribution(AminoAcid.PRO,AminoAcid.GLY);
        description = "";
        for (int i=0; i < dpdPro.outcomes.size(); i++)
            {
                Double probability = dpdPro.inputProbabilities.get(i) * 100.0;
                SidechainRotamerLibrary.BackboneAngles angles = dpdPro.outcomes.get(i);
                double phi = -1.0 * angles.getPhi();
                double psi = -1.0 * angles.getPsi();
                description += String.format("%6.1f,%6.1f,%9.6f\n", phi, psi, probability);
                if ( phi >= 51.0 && phi <= 108.0 && psi >= -10.0 && psi <= 47.0 )
                sum += probability;
            }
        System.out.println(sum);
        InputFileFormat.writeStringToDisk(description,"pro.csv");

        

        //System.out.println(getTripletDistribution(AminoAcid.LEU, AminoAcid.GLY, AminoAcid.MET).toDebugString(0.0001,10));
        
        //System.out.println(getTripletDistribution(AminoAcid.TYR, AminoAcid.ARG, AminoAcid.TRP).toDebugString(0.0001,10));
        //System.out.println(getTripletDistribution(AminoAcid.PHE, AminoAcid.ARG, AminoAcid.TRP).toDebugString(0.0001,10));
        
        //System.out.println(getTripletDistribution(AminoAcid.ASN, AminoAcid.TPR, AminoAcid.SER).toDebugString(0.0001,10));
        //System.out.println(getTripletDistribution(AminoAcid.ASN, AminoAcid.TPR, AminoAcid.THR).toDebugString(0.0001,10));
        Scanner scanner = new Scanner(System.in);
        System.out.println("Press enter to continue.");
        scanner.nextLine();
    }
}
