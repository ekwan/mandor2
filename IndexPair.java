import java.util.*;
import java.io.*;

/**
 * This is a utility generic class that represents a pair of rotamers by their indices.
 */
public class IndexPair implements Immutable, Serializable
{
    /** for serialization */
    public static final long serialVersionUID = 1L;

    /**
     * Represents the position on the odometer dial.
     * (i.e., the index of the residue in the sequence of the peptide.)
     */
    public final int outerIndex1, outerIndex2;

    /**
     * represents the number on the dial indicated by the outer index.
     * (i.e., which rotamer this is)
     */
    public final int innerIndex1, innerIndex2;

    /**
     * Constructs an immutable IndexPair
     */
    public IndexPair(int outerIndex1, int innerIndex1, int outerIndex2, int innerIndex2)
    {
        this.outerIndex1 = outerIndex1;
        this.innerIndex1 = innerIndex1;
        this.outerIndex2 = outerIndex2;
        this.innerIndex2 = innerIndex2;
    }

    /**
     * Reconstitutes the indexPair.
     * @param rotamerSpace the rotamer space the index pair refers to
     * @return the corresponding pair of rotamers
     */
    public RotamerPair getRotamerPair(List<List<Rotamer>> rotamerSpace)
    {
        Rotamer rotamer1 = rotamerSpace.get(outerIndex1).get(innerIndex1);
        Rotamer rotamer2 = rotamerSpace.get(outerIndex2).get(innerIndex2);
        return new RotamerPair(rotamer1, rotamer2);
    }

    /**
     * Determine if the given index is in this pair.  
     * @param outer The given outer index.  
     * @param inner The given inner index.  
     * @return true if either element of this pair equals the given index.  
     */
    public boolean hasPair(int outer, int inner)
    {
        return ((outerIndex1==outer&&innerIndex1==inner)||(outerIndex2==outer&&innerIndex2==inner));
    }

    public String toString()
    {
        return String.format("[%d, %d] : [%d, %d]", outerIndex1, innerIndex1, outerIndex2, innerIndex2);
    }

    /**
     * The minimum of the two outer indices that is not shared by given IndexPair.  
     * @param pair A pair of indices, the outer indices of which are not allowed as return values.  
     * @return The minimal outer index that isn't also an outer index of pair. 
     */
    public int minOuterIndex(IndexPair pair)
    {
        int lesser = ((outerIndex1<outerIndex2) ? outerIndex1 : outerIndex2);
        int greater = ((outerIndex1<outerIndex2) ? outerIndex2 : outerIndex1);
        if (lesser!=pair.outerIndex1&&lesser!=pair.outerIndex2) return lesser;
        return greater;
    }

    /**
     * Returns the hash code for this pair.
     * @return the hash code
     */
    @Override
    public int hashCode()
    {
        int result = 17;
        if ( outerIndex1 < outerIndex2 )
            {
                result = 31 * result + outerIndex1;
                result = 31 * result + innerIndex1;
                result = 31 * result + outerIndex2;
                result = 31 * result + innerIndex2;
            }
        else
            {
                result = 31 * result + outerIndex2;
                result = 31 * result + innerIndex2;
                result = 31 * result + outerIndex1;
                result = 31 * result + innerIndex1;
            }
        return result;
    }

    /**
     * Tests for object equality.
     * @return true if the pairs are equal
     */
    @Override
    public boolean equals(Object obj)
    {
        if ( obj == null )
            return false;
        if ( obj == this ) 
            return true;
        if ( !(obj instanceof IndexPair) )
            return false;

        IndexPair p = (IndexPair)obj;
        if ( outerIndex1 == p.outerIndex1 &&
             innerIndex1 == p.innerIndex1 &&
             outerIndex2 == p.outerIndex2 &&
             innerIndex2 == p.innerIndex2    )
            return true;
        else if ( outerIndex1 == p.outerIndex2 &&
                  innerIndex1 == p.innerIndex2 &&
                  outerIndex2 == p.outerIndex1 &&
                  innerIndex2 == p.innerIndex1    )
            return true;
        return false;
    }

    /** for testing */
    public static void main(String[] args)
    {
        IndexPair testPair = new IndexPair(1,2,3,4);
        System.out.println(testPair.hashCode());
        IndexPair testPair2 = new IndexPair(3,4,1,2);
        System.out.println(testPair2.hashCode());
        System.out.println(testPair.equals(testPair2));
    }
}
