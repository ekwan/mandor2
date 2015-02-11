import java.io.*;
import java.util.*;

/**
 * This represents a pair of rotamers.  Obeys (a,b).equals(b,a).
 */
public class RotamerPair implements Immutable
{
    /** one rotamer */
    public final Rotamer rotamer1;
    
    /** the other rotamer */
    public final Rotamer rotamer2;

    /**
     * Constructs a pair of rotamers.
     */
    public RotamerPair(Rotamer rotamer1, Rotamer rotamer2)
    {
        //if ( rotamer1 == null || rotamer2 == null )
        //    throw new NullPointerException("nulls not allowed");
        //if ( rotamer1.equals(rotamer2) )
        //    System.out.println("duplicate");
        this.rotamer1 = rotamer1;
        this.rotamer2 = rotamer2;
    }

    @Override
    public String toString()
    {
        return "=======RotamerPair========\n" + rotamer1.toString() + "\n" + rotamer2.toString() + "\n======";
    }

    /**
     * Returns the hash code for this pair.
     * @return the hash code
     */
    @Override
    public int hashCode()
    {
        if ( rotamer1.hashCode() < rotamer2.hashCode() )
            return Objects.hash(rotamer1, rotamer2);
        return Objects.hash(rotamer2, rotamer1);
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
        if ( !(obj instanceof RotamerPair) )
            return false;

        RotamerPair p = (RotamerPair)obj;
        if ( rotamer1.equals(p.rotamer1) && rotamer2.equals(p.rotamer2) )
            return true;
        else if ( rotamer1.equals(p.rotamer2) && rotamer2.equals(p.rotamer1) )
            return true;
        return false;
    }
}
