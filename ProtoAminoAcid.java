import java.util.*;
import java.io.*;
import com.google.common.collect.*;
import org.apache.commons.math3.geometry.euclidean.threed.*;
import java.util.concurrent.*;

/**
 * Represents a template geometry of an amino acid.
 * Extraneous atoms are not stripped.  This class is immutable.
 * Because of the way equals/hashCode works for atoms, it's not a good idea to have two different ProtoAminoAcids
 * that share common atoms.  Therefore, a shift(int) method is provided.
 */
public class ProtoAminoAcid implements Immutable, Serializable
{
    /** For serialization. */
    public static final long serialVersionUID = 1L;

    /** Defined as the ordered pair (terminal acetyl amide carbon, terminal acetyl amide nitrogen) */
    public final Pair<Atom,Atom> NStickyConnection;

    /** Defined as the ordered pair (amide carbonyl carbon, NH2 nitrogen) */
    public final Pair<Atom,Atom> CStickyConnection;

    /** Geometry of AcHN-Xxx-CONH2 exactly as read from the template file */
    public final Molecule molecule;

    /** Contains metadata about which atoms are which. */
    public final Residue residue;

    /** Constructs a ProtoAminoAcid from a file. */
    public ProtoAminoAcid(ProtoAminoAcidFile m)
    {
        this(m.NStickyConnection, m.CStickyConnection, m.molecule, m.residue);
   }

    /** Constructs a ProtoAminoAcid from a file. */
    private ProtoAminoAcid(Pair<Atom,Atom> NStickyConnection, Pair<Atom,Atom> CStickyConnection, Molecule molecule, Residue residue)
    {
        if ( NStickyConnection == null || CStickyConnection == null || molecule == null || residue == null )
            throw new NullPointerException("null fields are not allowed");
        
        this.molecule = molecule;
        this.NStickyConnection = NStickyConnection;
        this.CStickyConnection = CStickyConnection;
        if ( NStickyConnection.getFirst().element  != Element.NITROGEN ||
             NStickyConnection.getSecond().element != Element.CARBON   ||
             CStickyConnection.getFirst().element  != Element.CARBON   ||
             CStickyConnection.getSecond().element != Element.NITROGEN    )
            throw new IllegalArgumentException("unexpected sticky connection element");
        this.residue = residue;
        
        if ( ! molecule.contents.equals(residue.atoms) )
            throw new IllegalArgumentException("atom list mismatch");
    }

    /**
     * Returns a copy of this ProtoAminoAcid with the atoms shifted in coordinates by
     * (+i, +i, +i).
     * @param i how much to shift the geometry by (cannot be zero)
     * @return the shifted ProtoAminoAcid
     */
    public ProtoAminoAcid shift(int i)
    {
        if ( i == 0 )
            throw new IllegalArgumentException("nonzero shift required");
        Map<Atom,Atom> atomMap = new HashMap<>();
        for (Atom a : molecule.contents)
            {
                Vector3D oldPosition = a.position;
                Vector3D newPosition = oldPosition.add(new Vector3D(i*1.0, i*1.0, i*1.0));
                Atom newAtom = a.moveAtom(newPosition); 
                atomMap.put(a, newAtom);
            }
        return moveAtoms(atomMap);
    }

    /**
     * Creates a new ProtoAminoAcid from a map of old atoms to new atoms.
     * @param atomMap a map from old atoms to new atoms
     * @return the moved ProtoAminoAcid
     */
    public ProtoAminoAcid moveAtoms(Map<Atom,Atom> atomMap)
    {
        for (Atom a : atomMap.keySet())
            if ( ! molecule.contents.contains(a) )
                throw new IllegalArgumentException("atom map contains an atom that is not in this ProtoAminoAcid");

        Atom tempAtom1 = NStickyConnection.getFirst();
        Atom tempAtom2 = NStickyConnection.getSecond();
        if ( atomMap.containsKey(NStickyConnection.getFirst()) )
            tempAtom1 = atomMap.get(tempAtom1);
        if ( atomMap.containsKey(NStickyConnection.getSecond()) )
            tempAtom2 = atomMap.get(tempAtom2);
        Pair<Atom,Atom> newNStickyConnection = new Pair<>(tempAtom1,tempAtom2);

        tempAtom1 = CStickyConnection.getFirst();
        tempAtom2 = CStickyConnection.getSecond();
        if ( atomMap.containsKey(CStickyConnection.getFirst()) )
            tempAtom1 = atomMap.get(tempAtom1);
        if ( atomMap.containsKey(CStickyConnection.getSecond()) )
            tempAtom2 = atomMap.get(tempAtom2);
        Pair<Atom,Atom> newCStickyConnection = new Pair<>(tempAtom1,tempAtom2);

        Molecule newMolecule = molecule.moveAtoms(atomMap);
        Residue newResidue = residue.moveAtoms(atomMap);

        return new ProtoAminoAcid(newNStickyConnection, newCStickyConnection, newMolecule, newResidue);
    }

    @Override
    public String toString()
    {
        // print header, omega, phi, psi
        String returnString   = String.format("ProtoAminoAcid for %s\n\n", residue.aminoAcid.toString());
        returnString         += String.format("Description: %s\n", residue.description);
        returnString         += String.format("      Omega: %s\n", residue.omega.toString(molecule));
        returnString         += String.format("        Phi: %s\n", residue.phi.toString(molecule));
        returnString         += String.format("        Psi: %s\n", residue.psi.toString(molecule));

        // print chis
        for (int i=0; i < residue.chis.size(); i++)
            {
                ProtoTorsion t = residue.chis.get(i);
                returnString += String.format("     Chi %d: %s\n", i+1, t.toString(molecule));
            }
        if ( residue.chis.size() == 0 )
            returnString +=                   "       Chis: --none--\n";

        // print special atoms
        if ( residue.HN != null )
            returnString +=     String.format("         HN: %s\n", molecule.getAtomString(residue.HN));
        else
            returnString +=     String.format("         HN: none\n");
        returnString     +=     String.format("          N: %s\n", molecule.getAtomString(residue.N));
        returnString     +=     String.format("          O: %s\n", molecule.getAtomString(residue.O));
        returnString     +=     String.format("          C: %s\n", molecule.getAtomString(residue.C));
        returnString     +=     String.format("         CA: %s\n", molecule.getAtomString(residue.CA));
        if ( residue.HA != null )
            returnString +=     String.format("         HA: %s\n", molecule.getAtomString(residue.HA));
        else
            returnString +=     String.format("         HA: none\n");
 
        // print sticky connections
        returnString     += String.format("  NStickyConnection: %s-%s\n", molecule.getAtomString(NStickyConnection.getFirst()),
                                                                          molecule.getAtomString(NStickyConnection.getSecond()) );
        returnString     += String.format("  CStickyConnection: %s-%s\n", molecule.getAtomString(CStickyConnection.getFirst()),
                                                                          molecule.getAtomString(CStickyConnection.getSecond()) );
        returnString     += String.format("prochiralConnection: %s-%s",   molecule.getAtomString(residue.prochiralConnection.getFirst()),
                                                                          molecule.getAtomString(residue.prochiralConnection.getSecond()) );

        // return result
        return returnString;
    }

    @Override
    public int hashCode()
    {
        return Objects.hash(NStickyConnection, CStickyConnection, molecule, residue);
    }

    @Override
    public boolean equals(Object obj)
    {
        if ( obj == null )
            return false;
        if ( obj == this )
            return true;
        if ( !(obj instanceof ProtoAminoAcid) )
            return false;

        ProtoAminoAcid p = (ProtoAminoAcid)obj;
        if ( NStickyConnection.equals(p.NStickyConnection) &&
             CStickyConnection.equals(p.CStickyConnection) &&
             molecule.equals(p.molecule) &&
             residue.equals(p.residue) )
            return true;
        return false;
    }
}
