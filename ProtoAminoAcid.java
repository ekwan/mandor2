import java.util.*;
import java.io.*;
import com.google.common.collect.*;
import org.apache.commons.math3.geometry.euclidean.threed.*;
import java.util.concurrent.*;

/**
 * Represents a template geometry of an amino acid.
 * Extraneous atoms are not stripped.  This class is immutable.
 * Note that because I messed up when I designed everything, duplicate atoms in the template will
 * cause problems.  As a result, I provide methods for shifting the template atoms in space by
 * a random or deterministic amount.
 */
public class ProtoAminoAcid implements Immutable, Serializable
{
    /** for serialization */
    public static final long serialVersionUID = 1L;

    /** defined as the ordered pair (terminal acetyl amide carbon, terminal acetyl amide nitrogen) */
    public final Pair<Atom,Atom> NStickyConnection;

    /** defined as the ordered pair (amide carbonyl carbon, NH2 nitrogen) */
    public final Pair<Atom,Atom> CStickyConnection;

    /** geometry of AcHN-Xxx-CONH2 exactly as read from the template file */
    public final Molecule molecule;

    /** contains the residue information */
    public final Residue r;

    /** public constructor */
    public ProtoAminoAcid(MetadataOutputFile m)
    {
        Residue r = new Residue(m.aminoAcid, m.residueType, m.omega, m.phi, m.psi, m.chis, m.XHtorsion,
                                m.frozenTorsions, m.molecule.contents, m.frozenAtoms, m.backboneHN,
                                m.imidazoleHN, m.imidazoleN, m.otherHNatoms, m.referenceEnergy,
                                m.description, m.prochiralConnection, this);
        this.r = r;
        this.NStickyConnection = m.NStickyConnection;
        this.CStickyConnection = m.CStickyConnection;
        this.molecule = m.molecule;
    }

    private ProtoAminoAcid(Residue r, Molecule m, Pair<Atom,Atom> NStickyConnection, Pair<Atom,Atom> CStickyConnection)
    {
        this.r = r;
        this.molecule = m;
        this.NStickyConnection = NStickyConnection;
        this.CStickyConnection = CStickyConnection;
    }

    /**
     * Translates the atoms in this ProtoAminoAcid a random amount.
     * @return the shifted ProtoAminoAcid
     */
    public ProtoAminoAcid getRandomlyShifted()
    {
        // create random translation vector
        double deltaX = ThreadLocalRandom.current().nextDouble(-10.0, 10.0);
        double deltaY = ThreadLocalRandom.current().nextDouble(-10.0, 10.0);
        double deltaZ = ThreadLocalRandom.current().nextDouble(-10.0, 10.0);
        Vector3D translationVector = new Vector3D(deltaX, deltaY, deltaZ);

        // create atom map
        Map<Atom,Atom> atomMap = new HashMap<>();
        for (Atom oldAtom : molecule.contents)
            {
                Vector3D oldPosition = oldAtom.position;
                Vector3D newPosition = oldPosition.add(translationVector);
                Atom newAtom = oldAtom.moveAtom(newPosition);
                atomMap.put(oldAtom,newAtom);
            }

        // create new ProtoAminoAcid
        return moveAtoms(this, atomMap);
    }

    /**
     * Copy constructor that moves a ProtoAminoAcid.
     */
    public static ProtoAminoAcid moveAtoms(ProtoAminoAcid p, Map<Atom,Atom> atomMap)
    {
        Residue r  = p.r.moveAtoms(atomMap);
        Molecule m = p.molecule.moveAtoms(atomMap);

        Atom tempAtom1 = p.NStickyConnection.getFirst();
        Atom tempAtom2 = p.NStickyConnection.getSecond();
        if ( atomMap.containsKey(p.NStickyConnection.getFirst()) )
            tempAtom1 = atomMap.get(tempAtom1);
        if ( atomMap.containsKey(p.NStickyConnection.getSecond()) )
            tempAtom2 = atomMap.get(tempAtom2);
        Pair<Atom,Atom> NStickyConnection = new Pair<>(tempAtom1,tempAtom2);

        tempAtom1 = p.CStickyConnection.getFirst();
        tempAtom2 = p.CStickyConnection.getSecond();
        if ( atomMap.containsKey(p.CStickyConnection.getFirst()) )
            tempAtom1 = atomMap.get(tempAtom1);
        if ( atomMap.containsKey(p.CStickyConnection.getSecond()) )
            tempAtom2 = atomMap.get(tempAtom2);
        Pair<Atom,Atom> CStickyConnection = new Pair<>(tempAtom1,tempAtom2);

        return new ProtoAminoAcid(r, m, NStickyConnection, CStickyConnection);
    }

    /**
     * Returns a copy of this ProtoAminoAcid with the atoms shifted in coordinates by
     * (+i, +i, +i).
     */
    public ProtoAminoAcid shift(int i)
    {
        Map<Atom,Atom> atomMap = new HashMap<>();
        for (Atom a : molecule.contents)
            {
                Vector3D oldPosition = a.position;
                Vector3D newPosition = oldPosition.add(new Vector3D(i*1.0, i*1.0, i*1.0));
                Atom newAtom = a.moveAtom(newPosition); 
                atomMap.put(a, newAtom);
            }
        return moveAtoms(this,atomMap);
    }

    /**
     * Returns a textual description of all the ProtoAminoAcid fields.
     * @return the description
     */
    @Override
    public String toString()
    {
        // print header, omega, phi, psi
        String returnString =                String.format("ProtoAminoAcid for %s (%s)\n\n", r.aminoAcid.toString(), r.residueType.toString());
        returnString        = returnString + String.format("Description: %s\n", r.description);
        returnString        = returnString + String.format("Reference Energy: %.3f kcal/mol\n", r.referenceEnergy);
        returnString        = returnString + "Omega: " + r.omega.toString(molecule);
        if ( r.frozenTorsions.contains(r.omega) )
            returnString    = returnString + " (frozen)";
        returnString        = returnString + "\n";
        returnString        = returnString + "Phi:   " + r.phi.toString(molecule);
        if ( r.frozenTorsions.contains(r.phi) )
            returnString    = returnString + " (frozen)";
        returnString        = returnString + "\n";
        returnString        = returnString + "Psi  : " + r.psi.toString(molecule);
        if ( r.frozenTorsions.contains(r.psi) )
            returnString    = returnString + " (frozen)";
        returnString        = returnString + "\n";

        // print XH torsion
        if ( r.XHtorsion != null )
            returnString = returnString + "XHtorsion: " + r.XHtorsion.toString(molecule) + "\n"; 
        else
            returnString = returnString + "XHtorsion: none\n";
        
        // print chis
        for (int i=0; i < r.chis.size(); i++)
            {
                ProtoTorsion t = r.chis.get(i);
                returnString = returnString +    String.format("Chi %d: %s", i+1, t.toString(molecule));
                if ( r.frozenTorsions.contains(t) )
                    returnString = returnString + " (frozen)";
                returnString = returnString + "\n";
            }

        // print frozen atoms
        returnString = returnString + "Frozen atoms: ";
        for (Atom a : r.frozenAtoms)
            returnString = returnString + molecule.getAtomString(a) + ", ";
        if ( r.frozenAtoms.size() == 0 )
            returnString = returnString + "none";
        returnString = returnString + "\n";

        // print special atoms
        if ( r.backboneHN != null )
            returnString = returnString + String.format("Backbone HN: %s\n", molecule.getAtomString(r.backboneHN));
        else
            returnString = returnString + "Backbone HN: none\n";
        returnString = returnString + "Other HN atoms: ";
        for (Atom a : r.otherHNatoms)
            returnString = returnString + molecule.getAtomString(a) + ", ";
        if ( r.otherHNatoms.size() == 0 )
            returnString = returnString + "none";
        returnString = returnString + "\n";

        if ( r.aminoAcid == AminoAcid.HIS )
            {
                returnString = returnString + String.format("Imidazole HN: %s\n", molecule.getAtomString(r.imidazoleHN));
                returnString = returnString + String.format("Imidazole N:  %s\n", molecule.getAtomString(r.imidazoleN));
            }


        // print sticky connections
        returnString        = returnString + String.format("NStickyConnection:   %s-%s\n", molecule.getAtomString(NStickyConnection.getFirst()),
                                                                                           molecule.getAtomString(NStickyConnection.getSecond()) );
        returnString        = returnString + String.format("CStickyConnection:   %s-%s\n", molecule.getAtomString(CStickyConnection.getFirst()),
                                                                                           molecule.getAtomString(CStickyConnection.getSecond()) );
        
        returnString        = returnString + String.format("prochiralConnection: %s-%s",   molecule.getAtomString(r.prochiralConnection.getFirst()),
                                                                                           molecule.getAtomString(r.prochiralConnection.getSecond()) );

        // return result
        return returnString;
    }

    /**
     * Returns the hash code.
     * @return the hash code
     */
    @Override
    public int hashCode()
    {
        return Objects.hash(r, molecule, NStickyConnection, CStickyConnection);
    }

    /**
     * Tests object equality.
     * @param obj another object
     * @return true if equal
     */
    @Override
    public boolean equals(Object obj)
    {
        if ( obj == null )
            return false;
        if ( obj == this )
            return true;
        if ( !(obj instanceof Residue) )
            return false;

        ProtoAminoAcid p = (ProtoAminoAcid)obj;
        if ( Objects.equals(r, p.r) &&
             Objects.equals(molecule, p.molecule) && 
             Objects.equals(NStickyConnection, p.NStickyConnection) &&
             Objects.equals(CStickyConnection, p.CStickyConnection)    )
            return true;
        return false;
    }
}
