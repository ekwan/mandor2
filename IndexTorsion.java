import org.apache.commons.math3.geometry.euclidean.threed.*;
import java.io.*;
import java.util.*;
import com.google.common.collect.*;

/**
 * Represents a torsion composed of atom numbers in a molecule.<p>
 * This allows the IndexTorsion to remain usable if the molecular geometry changes. <p>
 * Torsion is: index1-index2-index3-index4<p>
 * Torsion angle calculated between atoms 2 and 3.<p>
 * If moved, the convention is to hold atom1 fixed and allow atom4 to move.<p>
 * This class contains a pre-computed list of the atoms to be rotated.
 */
public class IndexTorsion extends AbstractTorsion
{
    /** For serialization. */
    public static final long serialVersionUID = 1L;

    /** The atom indices (0,1,...,n-1) that comprise this torsion. */
    public final int index1, index2, index3, index4;

    /** The atom indices to rotate.  By convention, this is index4 and everything on that half of the connectivity graph. */
    public final List<Integer> atomIndicesToRotate;

    /** Constructs an IndexTorsion. */
    public IndexTorsion(int index1, int index2, int index3, int index4, List<Integer> atomIndicesToRotate)
    {
        this.index1 = index1;
        this.index2 = index2;
        this.index3 = index3;
        this.index4 = index4;
        
        // not checked for duplicates
        this.atomIndicesToRotate = ImmutableList.copyOf(atomIndicesToRotate);
        
        if ( index1 < 0 || index2 < 0 || index3 < 0 || index4 < 0 )
            throw new IllegalArgumentException("atom indices must be greater than zero");
        if ( ImmutableSet.of(index1,index2,index3,index4).size() != 4 )
            throw new IllegalArgumentException("duplicate atom indices");
        if ( atomIndicesToRotate == null || atomIndicesToRotate.size() == 0 )
            throw new NullPointerException("must have some atoms to rotate");
        if ( !atomIndicesToRotate.contains(index4) )
            throw new IllegalArgumentException("must rotate atom4");
        for (Integer i : atomIndicesToRotate)
            if ( i < 0 )
                throw new IllegalArgumentException("rotate atom indices must be greater than or equal to zero");
    }

    @Override
    public String toString()
    {
        return String.format("%d-%d-%d-%d: %s", index1, index2, index3, index4, atomIndicesToRotate.toString());
    }

    /**
     * Static factory method.
     */
    public static IndexTorsion createIndexTorsion(int index1, int index2, int index3, int index4, Molecule molecule)
    {
        Atom atom1 = molecule.contents.get(index1);
        Atom atom2 = molecule.contents.get(index2);
        Atom atom3 = molecule.contents.get(index3);
        Atom atom4 = molecule.contents.get(index4);
        if ( !molecule.directlyConnected(atom1,atom2) )
            throw new IllegalArgumentException("atom1-atom2 not connected");
        if ( !molecule.directlyConnected(atom2,atom3) )
            throw new IllegalArgumentException("atom1-atom2 not connected");
        if ( !molecule.directlyConnected(atom3,atom4) )
            throw new IllegalArgumentException("atom1-atom2 not connected");
        Set<Integer> atomIndicesToRotate = molecule.getHalfGraphIndices(index2, index3);
        return new IndexTorsion(index1,index2,index3,index4,ImmutableList.copyOf(atomIndicesToRotate));
    }

    /**
     * Makes an IndexTorsion from a ProtoTorsion.
     * @param torsion the ProtoTorsion to construct this IndexTorsion from
     * @param m the molecule the atoms in the ProtoTorsion belong to
     * @return the corresponding IndexTorsion
     */
    public static IndexTorsion createIndexTorsion(ProtoTorsion torsion, Molecule m)
    {
        int index1 = m.contents.indexOf(torsion.atom1);
        int index2 = m.contents.indexOf(torsion.atom2);
        int index3 = m.contents.indexOf(torsion.atom3);
        int index4 = m.contents.indexOf(torsion.atom4);

        if ( index1 == -1 || index2 == -1 || index3 == -1 || index4 == -1)
            throw new IllegalArgumentException("this atom is not in the molecule");

        return createIndexTorsion(index1,index2,index3,index4,m);
    }

    /** returns the corresponding torsion angle */
    public double getDihedralAngle()
    {
        throw new IllegalArgumentException("Mandor needs to know which Molecule this torsion corresponds to to calculate an angle.");
    }

    /**
     * Returns the corresponding torsion angle.  Could throw a NullPointerException.
     */
    public double getDihedralAngle(Molecule molecule)
    {
        Vector3D v1 = molecule.contents.get(index1).position; 
        Vector3D v2 = molecule.contents.get(index2).position; 
        Vector3D v3 = molecule.contents.get(index3).position; 
        Vector3D v4 = molecule.contents.get(index4).position;
        return getDihedralAngle(v1,v2,v3,v4);
    }

    /**
     * Returns the actual atoms to move in the specified molecule.
     * No checks.
     * @param molecule get the actual atoms from this molecule
     * @return the atoms to move on the atom4 side of the graph
     */
    public List<Atom> getAtomsToMove(Molecule molecule)
    {
        List<Atom> returnList = new ArrayList<>();
        for (Integer i : atomIndicesToRotate)
            returnList.add(molecule.contents.get(i));
        return ImmutableList.copyOf(returnList);
    }

    /**
     * Returns the corresponding ProtoTorsion.
     * @param molecule the molecule to resolve the indices on
     * @return the ProtoTorsion
     */
    public ProtoTorsion getProtoTorsion(Molecule molecule)
    {
        Atom atom1 = molecule.contents.get(index1);
        Atom atom2 = molecule.contents.get(index2);
        Atom atom3 = molecule.contents.get(index3);
        Atom atom4 = molecule.contents.get(index4);
        return new ProtoTorsion(atom1, atom2, atom3, atom4);
    }

    @Override
    public int hashCode()
    {
        return Objects.hash(index1, index2, index3, index4, atomIndicesToRotate);
    }

    @Override
    public boolean equals(Object obj)
    {
        if ( obj == null )
            return false;
        if ( obj == this )
            return true;
        if ( !(obj instanceof IndexTorsion) )
            return false;

        IndexTorsion t = (IndexTorsion)obj;
        if ( index1 == t.index1 &&
             index2 == t.index2 &&
             index3 == t.index3 &&
             index4 == t.index4 &&
             Objects.equals(atomIndicesToRotate, t.atomIndicesToRotate) )
            return true;
        return false;
    }
}
