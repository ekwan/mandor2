import org.apache.commons.math3.geometry.euclidean.threed.*;
import java.io.*;
import java.util.*;
import com.google.common.collect.*;

/**
 * Represents a torsion composed of concrete atoms.<p>
 * Torsion is: atom1-atom2-atom3-atom4<p>
 * Torsion angle calculated between atoms 2 and 3.<p>
 */
public class ProtoTorsion extends AbstractTorsion
{
    /** For serialization. */
    public static final long serialVersionUID = 1L; 

    /** The Atoms that comprise this torsion. */
    public final Atom atom1, atom2, atom3, atom4;

    /** Constructs a ProtoTorsion from four Atoms. */
    public ProtoTorsion(Atom atom1, Atom atom2, Atom atom3, Atom atom4)
    {
        if ( atom1 == null || atom2 == null || atom3 == null || atom4 == null )
            throw new NullPointerException("no nulls allowed");
        this.atom1 = atom1;
        this.atom2 = atom2;
        this.atom3 = atom3;
        this.atom4 = atom4;
    }

    /** Constructs a ProtoTorsion from four indices (0, 1, ..., n-1). */
    public ProtoTorsion(int index1, int index2, int index3, int index4, Molecule molecule)
    {
        this.atom1 = molecule.contents.get(index1);
        this.atom2 = molecule.contents.get(index2);
        this.atom3 = molecule.contents.get(index3);
        this.atom4 = molecule.contents.get(index4);
    }

    /** create a new ProtoTorsion using a map from old atoms to new atoms */
    public ProtoTorsion moveAtoms(Map<Atom,Atom> atomMap)
    {
        Atom newAtom1 = atom1;
        if ( atomMap.containsKey(atom1) )
            {
                newAtom1 = atomMap.get(atom1);
                //System.out.println(atom1 + " 1--> " + newAtom1);
            }
        
        Atom newAtom2 = atom2;
        if ( atomMap.containsKey(atom2) )
            {
                newAtom2 = atomMap.get(atom2);
                //System.out.println(atom2 + " 2--> " + newAtom2);
            }

        Atom newAtom3 = atom3;
        if ( atomMap.containsKey(atom3) )
            {
                newAtom3 = atomMap.get(atom3);
                //System.out.println(atom3 + " 3--> " + newAtom3);
            }

        Atom newAtom4 = atom4;
        if ( atomMap.containsKey(atom4) )
            {
                newAtom4 = atomMap.get(atom4);
                //System.out.println(atom4 + " 4--> " + newAtom4);
            }

        ProtoTorsion newProtoTorsion = new ProtoTorsion(newAtom1, newAtom2, newAtom3, newAtom4);
        return newProtoTorsion;
    }

    /** returns the corresponding torsion angle */
    public double getDihedralAngle()
    {
        return getDihedralAngle(atom1.position, atom2.position, atom3.position, atom4.position);
    }

    /**
     * Returns a string of the form C1-N2-H3-O4 given a molecule.
     * @param molecule the Molecule this ProtoTorsion corresponds to
     * @return the description
     */
    public String toString(Molecule molecule)
    {
        String atomString1 = molecule.getAtomString(atom1);
        String atomString2 = molecule.getAtomString(atom2);
        String atomString3 = molecule.getAtomString(atom3);
        String atomString4 = molecule.getAtomString(atom4);
        return String.format("%s-%s-%s-%s", atomString1, atomString2, atomString3, atomString4);
    }

    @Override
    public String toString()
    {
        return String.format("%s\n%s\n%s\n%s", atom1.position.toString(), atom2.position.toString(),
                                               atom3.position.toString(), atom4.position.toString()  );
    }

    /**
     * Returns the hash code for this AtomTorsion
     * @return the hash code
     */
    @Override
    public int hashCode()
    {
        return Objects.hash(atom1, atom2, atom3, atom4);
    }

    /**
     * Returns true if the AtomTorsions are equal.
     * @return true if equal
     */
    @Override
    public boolean equals(Object obj)
    {
        if ( obj == null )
            return false;
        if ( obj == this )
            return true;
        if ( !(obj instanceof ProtoTorsion) )
            return false;

        ProtoTorsion t = (ProtoTorsion)obj;
        if ( Objects.equals(atom1, t.atom1) &&
             Objects.equals(atom2, t.atom2) &&
             Objects.equals(atom3, t.atom3) &&
             Objects.equals(atom4, t.atom4)    )
            return true;
        return false;
    }
}
