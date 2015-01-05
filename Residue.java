import java.util.*;
import java.io.*;
import com.google.common.collect.*;

/**
 * This class represents an amino acid within a Peptide.
 * This class is immutable.
 */
public class Residue implements Immutable, Serializable
{
    /** for serialization */
    public static final long serialVersionUID = 1L;

    /** the kind of amino acid this is */
    public final AminoAcid aminoAcid;

    /** chirality if this is a normal amino acid or if this is a special transition state amino acid */
    public final ResidueType residueType;

    /** omega torsion */
    public final ProtoTorsion omega;

    /** phi torsion */
    public final ProtoTorsion phi;

    /** psi torsion */
    public final ProtoTorsion psi;
    
    /** chi1, chi2, chi3, chi4 */
    public final List<ProtoTorsion> chis;

    /** the atoms that are in this residue for energy calculations */
    public final List<Atom> atoms;

    /** the backbone amide NH */
    public final Atom backboneHN;

    /** if this is an imidazole, the imidazole NH; null otherwise */
    public final Atom imidazoleHN;

    /** if this is an imidazole, the imidazole N; null otherwise */
    public final Atom imidazoleN;

    /** any other amide NHs */
    public final List<Atom> otherHNatoms;

    /** reference energy in kcal/mol for comparing peptides of different composition */
    public final double referenceEnergy;

    /** the description of this residue */
    public final String description;

    /** the bond between the alpha carbon and the pro-L atom (usually the beta carbon) */
    public final Pair<Atom,Atom> prochiralConnection;

    /** public constructor */
    public Residue(AminoAcid aminoAcid, ResidueType residueType, ProtoTorsion omega,
                   ProtoTorsion phi, ProtoTorsion psi, List<ProtoTorsion> chis, List<Atom> atoms,
                   Atom backboneHN, Atom imidazoleHN,
                   Atom imidazoleN, List<Atom> otherHNatoms, double referenceEnergy,
                   String description, Pair<Atom,Atom> prochiralConnection)
    {
        // check invariants
        if ( aminoAcid == null || residueType == null || omega == null || phi == null ||
             psi == null || chis == null || frozenTorsions == null || atoms == null ||
             frozenAtoms == null || otherHNatoms == null )
            throw new NullPointerException("nulls not allowed");
        if ( atoms.size() == 0 )
            throw new IllegalArgumentException("atom list is empty");

        this.aminoAcid = aminoAcid;
        this.residueType = residueType;
        this.omega = omega;
        this.phi = phi;
        this.psi = psi;
        this.chis = ImmutableList.copyOf(chis);
        this.atoms = ImmutableList.copyOf(atoms);
        this.backboneHN = backboneHN;
        this.imidazoleHN = imidazoleHN;
        this.imidazoleN = imidazoleN;
        this.otherHNatoms = ImmutableList.copyOf(otherHNatoms);
        this.referenceEnergy = referenceEnergy;
        this.description = description;
        this.prochiralConnection = prochiralConnection;
    }

    
    /** returns a new residue with an updated atom list. Meant for adding HN if switching from proline in SidechainMutator */
    public Residue getNewResidue(List<Atom> newAtoms, Atom backboneHN)
    {
        //Add backbone HN
        return new Residue(aminoAcid, residueType, omega, phi, psi, chis, newAtoms, backboneHN, imidazoleHN, imidazoleN, otherHNatoms, referenceEnergy, description,prochiralConnection);
    }

    /** returns a new residue with an updated atom list. Meant for removing HN if switching to proline in SidechainMutator */
    public Residue getNewResidue(List<Atom> newAtoms)
    {
        //Remove backbone HN
        return new Residue(aminoAcid, residueType, omega, phi, psi, chis, newAtoms, null, imidazoleHN, imidazoleN, otherHNatoms, referenceEnergy, description,prochiralConnection);
    }

    /** returns a new residue given a map of old atoms to new atoms */
    public Residue moveAtoms(Map<Atom,Atom> atomMap)
    {
        ProtoTorsion newOmega                    = omega.moveAtoms(atomMap);
        ProtoTorsion newPhi                      = phi.moveAtoms(atomMap);
        ProtoTorsion newPsi                      = psi.moveAtoms(atomMap);
        
        List<ProtoTorsion> newChis               = new LinkedList<>();
        for (ProtoTorsion t : chis)
            newChis.add( t.moveAtoms(atomMap) );
        newChis = ImmutableList.copyOf(newChis);
      
        List<Atom> newAtoms                      = copyList(atoms, atomMap);

        Atom newBackboneHN                       = null;
        if ( backboneHN != null )
            newBackboneHN = backboneHN.moveAtom(atomMap);

        Atom newImidazoleHN                      = null;
        Atom newImidazoleN                       = null;
        if ( aminoAcid == AminoAcid.HIS )
            {
                newImidazoleHN = imidazoleHN.moveAtom(atomMap);
                newImidazoleN  = imidazoleN.moveAtom(atomMap);
            }

        List<Atom> newOtherHNatoms               = copyList(otherHNatoms, atomMap);

        Atom tempAtom1 = prochiralConnection.getFirst();
        Atom tempAtom2 = prochiralConnection.getSecond();
        if ( atomMap.containsKey(tempAtom1) )
            tempAtom1 = atomMap.get(tempAtom1);
        if ( atomMap.containsKey(tempAtom2) )
            tempAtom2 = atomMap.get(tempAtom2);
        Pair<Atom,Atom> newProchiralConnection   = new Pair<Atom,Atom>(tempAtom1, tempAtom2);

        return new Residue(aminoAcid, residueType, newOmega, newPhi, newPsi,
                           newChis, newAtoms, newBackboneHN, newImidazoleHN, newImidazoleN,
                           newOtherHNatoms, referenceEnergy, description, newProchiralConnection);
    }

    /** copies a list of atoms to a new list of atoms given a map from old atoms to old */
    private static List<Atom> copyList(List<Atom> list, Map<Atom,Atom> atomMap)
    {
        List<Atom> newList = new LinkedList<>();
        for (Atom a : list)
            {
                if ( atomMap.containsKey(a) )
                    newList.add(atomMap.get(a));
                else
                    newList.add(a);
            }
        return ImmutableList.copyOf(newList);
    }

    /**
     * Returns a description of this Residue.
     * @return the description
     */
    @Override
    public String toString()
    {
        return String.format("Residue %s (%s, %s)", aminoAcid.toString(), residueType.toString(), description);
    }

    /**
     * Returns a detailed description of this Residue.
     * @param peptide the peptide this residue is associated with
     * @return the textual description
     */
    public String toString(Peptide peptide)
    {
        String returnString = String.format("%s, %s, %s (%d atoms, ref energy = %.1f)\n", aminoAcid.toString(), residueType.toString(),
                                            description, atoms.size(), referenceEnergy);
        returnString = returnString + String.format(" omega: %-20s %7.1f\n", omega.toString(peptide), omega.getDihedralAngle());
        returnString = returnString + String.format("   phi: %-20s %7.1f\n", phi.toString(peptide), phi.getDihedralAngle());
        returnString = returnString + String.format("   psi: %-20s %7.1f\n", psi.toString(peptide), psi.getDihedralAngle());

        for (ProtoTorsion p : chis)
            returnString = returnString + String.format("   chi: %-20s %7.1f\n", p.toString(peptide), p.getDihedralAngle());
        
        if ( backboneHN != null )
            returnString = returnString + "backboneHN: " + peptide.getAtomString(backboneHN) + "\n";
        else
            returnString = returnString + "backboneHN: null\n";

        if ( imidazoleHN != null )
            returnString = returnString + "imidazoleHN: " + peptide.getAtomString(imidazoleHN) + "\n";
        else
            returnString = returnString + "imidazoleHN: null\n";

        if ( imidazoleN != null )
            returnString = returnString + "imidazoleN: " + peptide.getAtomString(imidazoleN) + "\n";
        else
            returnString = returnString + "imidazoleN: null\n";

        returnString = returnString + "Other HNatoms: ";
        for ( Atom a : otherHNatoms )
            returnString = returnString + peptide.getAtomString(a) + " ";
        
        returnString = returnString + "\n";

        returnString = returnString + String.format("prochiral connection: %s-%s\n", peptide.getAtomString(prochiralConnection.getFirst()),
                                                    peptide.getAtomString(prochiralConnection.getSecond()) );
        return returnString;
    }

    /**
     * Returns the hash code.  To prevent circular hash codes, protoAminoAcid is not considered.
     * @return the hash code
     */
    @Override
    public int hashCode()
    {
        return Objects.hash(aminoAcid, residueType, omega, phi, psi, chis, atoms, backboneHN,
                            imidazoleHN, imidazoleN, otherHNatoms, referenceEnergy,
                            description, prochiralConnection);
    }

    /**
     * Tests for object equality.  To prevent circular hash codes, protoAminoAcid is not considered.
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

        Residue r = (Residue)obj;
        if ( aminoAcid == r.aminoAcid &&
             residueType == r.residueType &&
             Objects.equals(omega, r.omega) &&
             Objects.equals(phi, r.phi) &&
             Objects.equals(psi, r.psi) &&
             Objects.equals(chis, r.chis) &&
             Objects.equals(atoms, r.atoms) &&
             Objects.equals(backboneHN, r.backboneHN) &&
             Objects.equals(imidazoleHN, r.imidazoleHN) &&
             Objects.equals(imidazoleN, r.imidazoleN) &&
             Objects.equals(otherHNatoms, r.otherHNatoms) &&
             (referenceEnergy == r.referenceEnergy) &&
             Objects.equals(description, r.description) &&
             Objects.equals(prochiralConnection, r.prochiralConnection) )
            return true;
        return false;
    }
}
