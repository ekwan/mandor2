import java.util.*;
import java.io.*;
import com.google.common.collect.*;

/**
 * This immutable class represents an amino acid within a Peptide.
 *
 * A residue is this chemical fragment:
 *
 * --HN-CO-CH(sidechain)--
 *
 * Fields for the atoms are stored using PDB nomenclature.
 * All fields are mandatory unless otherwise specified.
 */
public class Residue implements Immutable, Serializable
{
    /** for serialization */
    public static final long serialVersionUID = 1L;

    /** the kind of amino acid this is */
    public final AminoAcid aminoAcid;

    /** omega torsion */
    public final ProtoTorsion omega;

    /** phi torsion */
    public final ProtoTorsion phi;

    /** psi torsion */
    public final ProtoTorsion psi;
    
    /** chi torsions (if no torsions, use an empty list) */
    public final List<ProtoTorsion> chis;

    /** the backbone amide hydrogen (optional) */
    public final Atom HN;

    /** the backbone amide nitrogen */
    public final Atom N;

    /** the carbonyl oxygen */
    public final Atom O;

    /** the carbonyl carbon */
    public final Atom C;

    /** the alpha carbon */
    public final Atom CA;

    /** the alpha hydrogen (optional) */
    public final Atom HA;

    /** the description of this residue */
    public final String description;

    /** the bond between the alpha carbon and the pro-L atom (usually the beta carbon) */
    public final Pair<Atom,Atom> prochiralConnection;

    /** the atoms in this residue */
    public final List<Atom> atoms;

    /** whether this is a hairpin residue or not -- note this is initially set to false, use PeptideFactory.setHairpinAngles */
    public final boolean isHairpin;

    /** Constructs a Residue.  Many things are checked, but the checks are not exhaustive. */
    public Residue(AminoAcid aminoAcid, ProtoTorsion omega, ProtoTorsion phi, ProtoTorsion psi, List<ProtoTorsion> chis,
                   Atom HN, Atom N, Atom O, Atom C, Atom CA, Atom HA, String description,
                   Pair<Atom,Atom> prochiralConnection, List<Atom> atoms, boolean isHairpin) 
    {
        if ( aminoAcid == null )
            throw new NullPointerException("null amino acid");
        this.aminoAcid = aminoAcid;

        if ( omega == null )
            throw new NullPointerException("null omega");
        this.omega = omega;

        if ( phi == null )
            throw new NullPointerException("null phi");
        this.phi = phi;

        if ( psi == null )
            throw new NullPointerException("null psi");
        this.psi = psi;

        if ( chis == null )
            throw new NullPointerException("null chis");
        this.chis = ImmutableList.copyOf(chis);

        this.HN = HN;
        if ( HN != null && HN.element != Element.HYDROGEN )
            throw new IllegalArgumentException("unexpected element for HN");

        if ( N == null )
            throw new NullPointerException("null N");
        if ( N.element != Element.NITROGEN )
            throw new IllegalArgumentException("unexpected element for N");
        this.N = N;

        if ( O == null )
            throw new NullPointerException("null O");
        if ( O.element != Element.OXYGEN )
            throw new IllegalArgumentException("unexpected element for O");
        this.O = O;

        if ( C == null )
            throw new NullPointerException("null C");
        if ( C.element != Element.CARBON )
            throw new IllegalArgumentException("unexpected element for C");
        this.C = C;

        if ( CA == null )
            throw new NullPointerException("null CA");
        if ( CA.element != Element.CARBON )
            throw new IllegalArgumentException("unexpected element for CA");
        this.CA = CA;

        this.HA = HA;
        if ( HA != null && HA.element != Element.HYDROGEN )
            throw new IllegalArgumentException("unexpected element for HA");

        if ( description == null || description.trim().length() == 0 )
            throw new NullPointerException("null or zero length description");
        this.description = description;

        if ( prochiralConnection == null )
            throw new NullPointerException("null prochiral connection");
        if ( !prochiralConnection.getFirst().equals(CA) )
            throw new IllegalArgumentException("first atom of prochiral connection must be the CA");
        this.prochiralConnection = prochiralConnection;

        if ( atoms == null )
            throw new NullPointerException("null atoms");
        if ( atoms.size() == 0 )
            throw new IllegalArgumentException("zero size atoms");
        
        // check that all the atoms in the fields are contained in the atom list
        // note that the residue atoms might not contain all of the atoms in the prototorsions because they might extend out of the residue
        Set<Atom> tempSet = ImmutableSet.of(N, O, C, CA, prochiralConnection.getFirst(), prochiralConnection.getSecond(),
                                            omega.atom3, omega.atom4, phi.atom2, phi.atom3, phi.atom4,
                                            psi.atom1, psi.atom2, psi.atom3);
        tempSet = new TreeSet<>(tempSet);
        for (ProtoTorsion t : chis)
            tempSet.addAll(ImmutableSet.of(t.atom1, t.atom2, t.atom3, t.atom4));
        if ( HN != null )
            tempSet.add(HN);
        if ( HA != null )
            tempSet.add(HA);
        if ( ! atoms.containsAll(tempSet) )
            {
                System.out.println("missing:");
                for (Atom a : tempSet)
                    {
                        if ( !atoms.contains(a) )
                            System.out.println(a.toFullString());
                    }
                throw new IllegalArgumentException("atom list does not contain all atoms in fields");
            }

        // check that the torsions line up with the fields
        if ( ImmutableSet.of(omega.atom2, phi.atom1).size() != 1 )
            throw new IllegalArgumentException("carbon prior to N doesn't line up");
        if ( ImmutableSet.of(N, omega.atom3, phi.atom2, psi.atom1).size() != 1 )
            throw new IllegalArgumentException("N does not line up");
        if ( ImmutableSet.of(CA, omega.atom4, phi.atom3, psi.atom2).size() != 1 )
            throw new IllegalArgumentException("CA does not line up");
        if ( ImmutableSet.of(C, phi.atom4, psi.atom3).size() != 1 )
            throw new IllegalArgumentException("C does not line up");
        
        if (chis.size() == 1)
            {
                ProtoTorsion chi1 = chis.get(0);
                if ( ImmutableSet.of(N,  chi1.atom1).size() != 1 ||
                     ImmutableSet.of(CA, chi1.atom2).size() != 1    )
                    throw new IllegalArgumentException("problem with chi1");
            }
        else if ( chis.size() == 2 )
            {
                ProtoTorsion chi1 = chis.get(0);
                ProtoTorsion chi2 = chis.get(1);
                if ( ImmutableSet.of(N,  chi1.atom1).size() != 1 ||
                     ImmutableSet.of(CA, chi1.atom2).size() != 1 )
                    throw new IllegalArgumentException("problem with chi1");
                if ( ImmutableSet.of(CA, chi2.atom1).size() != 1 )
                    throw new IllegalArgumentException("problem with chi2");
                if ( ImmutableSet.of(chi1.atom2, chi2.atom1).size() != 1 ||
                     ImmutableSet.of(chi1.atom3, chi2.atom2).size() != 1 ||
                     ImmutableSet.of(chi1.atom4, chi2.atom3).size() != 1    )
                    throw new IllegalArgumentException("chis do not line up");
            }
        else if ( chis.size() == 3 )
            {
                ProtoTorsion chi1 = chis.get(0);
                ProtoTorsion chi2 = chis.get(1);
                ProtoTorsion chi3 = chis.get(2);
                if ( ImmutableSet.of(N,  chi1.atom1).size() != 1 ||
                     ImmutableSet.of(CA, chi1.atom2).size() != 1 )
                    throw new IllegalArgumentException("problem with chi1");
                if ( ImmutableSet.of(CA, chi2.atom1).size() != 1 )
                    throw new IllegalArgumentException("problem with chi2");
                if ( ImmutableSet.of(chi1.atom2, chi2.atom1).size() != 1 ||
                     ImmutableSet.of(chi1.atom3, chi2.atom2, chi3.atom1).size() != 1 ||
                     ImmutableSet.of(chi1.atom4, chi2.atom3, chi3.atom2).size() != 1 ||
                     ImmutableSet.of(chi2.atom4, chi3.atom3).size() != 1                )
                    throw new IllegalArgumentException("chis do not line up");
            }
        else if ( chis.size() == 4 )
             {
                ProtoTorsion chi1 = chis.get(0);
                ProtoTorsion chi2 = chis.get(1);
                ProtoTorsion chi3 = chis.get(2);
                ProtoTorsion chi4 = chis.get(3);
                if ( ImmutableSet.of(N,  chi1.atom1).size() != 1 ||
                     ImmutableSet.of(CA, chi1.atom2).size() != 1 )
                    throw new IllegalArgumentException("problem with chi1");
                if ( ImmutableSet.of(CA, chi2.atom1).size() != 1 )
                    throw new IllegalArgumentException("problem with chi2");
                if ( ImmutableSet.of(chi1.atom2, chi2.atom1).size() != 1 ||
                     ImmutableSet.of(chi1.atom3, chi2.atom2, chi3.atom1).size() != 1 ||
                     ImmutableSet.of(chi1.atom4, chi2.atom3, chi3.atom2, chi4.atom1).size() != 1 ||
                     ImmutableSet.of(chi2.atom4, chi3.atom3, chi4.atom2).size() != 1 ||
                     ImmutableSet.of(chi3.atom4, chi4.atom3).size() != 1 )
                    throw new IllegalArgumentException("chis do not line up");
            }
        else if ( chis.size() > 4 )
            throw new IllegalArgumentException("don't know how to deal with more than four chis");

        this.atoms = ImmutableList.copyOf(atoms);
        this.isHairpin = isHairpin;
    }

    /** Returns a new residue given a map of old atoms to new atoms. */
    public Residue moveAtoms(Map<Atom,Atom> atomMap)
    {
        // move torsions
        ProtoTorsion newOmega = omega.moveAtoms(atomMap);
        ProtoTorsion newPhi   = phi.moveAtoms(atomMap);
        ProtoTorsion newPsi   = psi.moveAtoms(atomMap);
        List<ProtoTorsion> newChis = new ArrayList<>();
        for (ProtoTorsion t : chis)
            newChis.add(t.moveAtoms(atomMap));

        // move atom fields
        Atom newHN = null;
        if ( HN != null )
            newHN = HN.moveAtom(atomMap);
        Atom newN  = N.moveAtom(atomMap);
        Atom newO  = O.moveAtom(atomMap);
        Atom newC  = C.moveAtom(atomMap);
        Atom newCA = CA.moveAtom(atomMap);
        Atom newHA = null;
        if ( HA != null )
            newHA = HA.moveAtom(atomMap);

        // move the prochiral connection
        Atom tempAtom1 = prochiralConnection.getFirst();
        Atom tempAtom2 = prochiralConnection.getSecond();
        if ( atomMap.containsKey(tempAtom1) )
            tempAtom1 = atomMap.get(tempAtom1);
        if ( atomMap.containsKey(tempAtom2) )
            tempAtom2 = atomMap.get(tempAtom2);
        Pair<Atom,Atom> newProchiralConnection   = new Pair<Atom,Atom>(tempAtom1, tempAtom2);
    
        // move the atom list
        List<Atom> newAtoms = new ArrayList<>();
        for (Atom a : atoms)
            newAtoms.add(a.moveAtom(atomMap));

        // return result
        return new Residue(aminoAcid, newOmega, newPhi, newPsi, newChis, newHN, newN, newO, newC, newCA, newHA, description, newProchiralConnection, newAtoms, isHairpin);
    }

    /** 
     * Returns the same residue with a different hairpin flag.
     * @param newIsHairpin whether this is a hairpin residue
     * @return the new residue
     */
    public Residue setHairpin(boolean newIsHairpin)
    {
        if ( newIsHairpin == isHairpin)
            return this;
        return new Residue(aminoAcid, omega, phi, psi, chis, HN, N, O, C, CA, HA, description, prochiralConnection, atoms, newIsHairpin);
    }

    @Override
    public String toString()
    {
        return String.format("Residue %s (%s)", aminoAcid.toString(), description);
    }

    /**
     * Returns a detailed description of this Residue.
     * @param peptide the peptide this residue is associated with
     * @return the textual description
     */
    public String toString(Molecule peptide)
    {
        String returnString = String.format("%s (%s, %d atoms)\n", aminoAcid.toString(), description, atoms.size());
        returnString += String.format("      omega: %-20s %7.1f\n", omega.toString(peptide), omega.getDihedralAngle());
        returnString += String.format("        phi: %-20s %7.1f\n", phi.toString(peptide), phi.getDihedralAngle());
        returnString += String.format("        psi: %-20s %7.1f\n", psi.toString(peptide), psi.getDihedralAngle());

        for (ProtoTorsion p : chis)
            returnString += String.format("        chi: %-20s %7.1f\n", p.toString(peptide), p.getDihedralAngle());
        
        List<Atom> tempList = new ArrayList<>();
        tempList.add(HN);
        tempList.add(N);
        tempList.add(O);
        tempList.add(C);
        tempList.add(CA);
        tempList.add(HA);
        returnString += "   HN,    N,    O,    C,   CA,   HA\n";
        for (Atom a : tempList)
            {
                if ( a != null )
                    returnString += String.format("%5s,", peptide.getAtomString(a));
                else
                    returnString += "  ---,";
            } 
        returnString = returnString.substring(0, returnString.length()-1);
        returnString += "\n";

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
        return Objects.hash(aminoAcid, atoms);
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
             atoms.equals(r.atoms) )
            return true;
        return false;
    }
}
