import org.apache.commons.math3.geometry.euclidean.threed.*;
import java.io.*;
import java.util.*;
import com.google.common.collect.*;

/**
 *  Represents an atom.  This class is immutable.
 */
public class Atom implements Immutable, Serializable, Comparable<Atom>
{
    /** For serialization. */
    public static final long serialVersionUID = 1L;

    /** Atomic element. */
    public final Element element;

    /** Location of the atom. */
    public final Vector3D position;

    /** AMOEBA atom type. */
    public final int type1;

    /** OPLS atom type. */
    public final int type2;

    /** The surface tension of this atom in kcal/angstroms^2. */
    public final double surfaceTension;

    /**
     * Constructs a new atom.
     * @param symbol the atomic symbol
     * @param position the location of the atom
     * @param type1 the AMOEBA atom type
     * @param type2 the OPLS atom type
     * @param surfaceTension the surface tension of this atom in kcal/angstroms^2
     */
    public Atom(String symbol, Vector3D position, int type1, int type2, double surfaceTension)
    {
        this(Element.getElement(symbol), position, type1, type2, surfaceTension);
    }

    /** Constructs a new Atom. */
    public Atom(Element element, Vector3D position, int type1, int type2, double surfaceTension)
    {
        this.element  = element;
        this.position = position;
        if ( type1 < 0 || type2 < 0 )
            throw new IllegalArgumentException("negative atom type");
        this.type1    = type1;
        this.type2    = type2;
        this.surfaceTension = surfaceTension;
    }   

    /** 
     * Compares the location of one atom to another using a comparison chart.
     * @param a the atom being compared to this atom
     * @return an int representing the result of the comparison (1 - same, 0 - different location) 
     */
    @Override
    public int compareTo(Atom a)
    {
        double x1 =   position.getX();
        double x2 = a.position.getX();
        double y1 =   position.getY();
        double y2 = a.position.getY();
        double z1 =   position.getZ();
        double z2 = a.position.getZ();
        return ComparisonChain.start().compare(x1, x2).compare(y1, y2).compare(z1, z2).result();
    }

    /**
     * Returns a copy of this atom with a new atom type.
     * @param newType1 the AMOEBA atom type for the new atom
     * @param newType2 the OPLS atom type for the new atom
     * @return a new atom with updated atom types
     */
    public Atom setAtomType(int newType1, int newType2)
    {
        return new Atom(element, position, newType1, newType2, surfaceTension);
    }

    public Atom changeTypes(Atom anotherAtom)
    {
        return new Atom(element, position, anotherAtom.type1, anotherAtom.type2, anotherAtom.surfaceTension);
    }

    /**
     * Returns a copy of this atom with a new position.
     * @param newPosition the new position
     * @return the new atom
     */
    public Atom moveAtom(Vector3D newPosition)
    {
        return new Atom(element.symbol, newPosition, type1, type2, surfaceTension);
    }

    /**
     * Rotates and translates this atom.  The rotation is applied before the translation.  
     * @param rot a three-dimensional rotation.
     * @param shift a vector that we add to the position of this after it has been rotated.
     * @return the new atom
     */
    public Atom transform(Rotation rot, Vector3D shift)
    {
        return new Atom(element.symbol, rot.applyTo(position).add(shift), type1, type2, surfaceTension);
    }

    /**
     * Convenience method for moving an atom based on a map of old atoms and new atoms.
     * @param atomMap a map from old atoms to new atoms
     * @return the new atom if applicable
     */
    public Atom moveAtom(Map<Atom,Atom> atomMap)
    {
        if (atomMap.containsKey(this))
            return atomMap.get(this);
        else
            return this;
    }

    @Override
    public String toString()
    {
        return String.format("%-2s %10.6f %10.6f %10.6f", element.symbol, position.getX(), position.getY(), position.getZ());
    }

    public String toFullString()
    {
        return String.format("%-2s %10.6f %10.6f %10.6f %3d %3d %7.4f", element.symbol, position.getX(), position.getY(), position.getZ(), type1, type2, surfaceTension);
    }


    @Override
    public int hashCode()
    {
        return Objects.hash(element, position, type1, type2, surfaceTension);
    }

    @Override
    public boolean equals(Object obj)
    {
        if ( obj == null )
            return false;
        else if ( obj == this )
            return true;
        else if ( !(obj instanceof Atom) )
            return false;

        Atom a = (Atom)obj;
        if ( element == a.element &&
             position.equals(a.position) &&
             type1 == a.type1 &&
             type2 == a.type2 &&
             surfaceTension == a.surfaceTension )
            return true;
        return false;
    }

    /** For testing. */
    public static void main(String[] args)
    {
        Atom atom1 = new Atom("C", new Vector3D(1.0, 2.0, 3.0), 103, 10, 0.0);
        Atom atom2 = new Atom("C", new Vector3D(2.0, 2.0, 3.0), 103, 10, 0.0);
        Atom atom3 = new Atom("C", new Vector3D(3.0, 3.0, 3.0), 103, 10, 0.0);
        Atom atom4 = new Atom("C", new Vector3D(3.0, 2.0, 3.0), 103, 10, 0.0);
        Atom atom5 = new Atom("C", new Vector3D(3.0, 2.0, 4.0), 103, 10, 0.0);
        List<Atom> list = new LinkedList<>();
        list.add(atom2);
        list.add(atom1);
        list.add(atom3);
        list.add(atom4);
        list.add(atom5);
        Collections.shuffle(list);
        Collections.sort(list);
        for (Atom a : list)
            System.out.println(a);
    }
} // end of class Atom
