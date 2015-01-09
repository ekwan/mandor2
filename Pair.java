import java.io.*;
import java.util.*;

/**
 * This is a utility generic class that represents an ordered pair.
 * This class is immutable.
 */
public class Pair<K,V> implements Serializable, Immutable
{
    public static final long serialVersionUID = 1L;

    private final K firstValue;
    private final V secondValue;
    
    /**
     * Constructs an immutable pair.  Nulls are not allowed.
     */
    public Pair(K firstVal, V secondVal)
    {
        if ( firstVal == null || secondVal == null )
            throw new NullPointerException("nulls not allowed");
	    this.firstValue = firstVal;
	    this.secondValue = secondVal;
    }

    /**
     * Returns the first value of this pair.
     * @return the first value
     */
    public K getFirst()
    {
	    return firstValue;
    }

    /**
     * Returns the second value of this pair.
     * @return the second value
     */
    public V getSecond()
    {
	    return secondValue;
    }

    /**
     * Returns a short description of this pair.
     * @return [first value, second value]
     */
    public String toString()
    {
        return String.format("[%s, %s]", firstValue.toString(), secondValue.toString());
    }

    /**
     * Returns the hash code for this pair.
     * @return the hash code
     */
    @Override
    public int hashCode()
    {
        return Objects.hash(firstValue, secondValue);
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
        if ( !(obj instanceof Pair) )
            return false;

        Pair<?,?> anotherPair = (Pair<?,?>)obj;
        if ( this.firstValue.equals(anotherPair.firstValue) &&
             this.secondValue.equals(anotherPair.secondValue)  )
            return true;
        return false;
    }
}
