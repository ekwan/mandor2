import org.apache.commons.math3.geometry.euclidean.threed.*;
import java.io.*;
import java.util.*;
import java.util.concurrent.*;
import org.jgrapht.*;
import org.jgrapht.graph.*;
import org.jgrapht.alg.*;
import com.google.common.collect.*;

/**
 * This collects together methods for taking a peptide and changing its
 * the sidechain chi angles.  Dihedral angles are chosen in a backbone-dependent
 * way; i.e., the residue's (phi,psi) determines the conditional probability
 * disttribution for the {chi_N}.
 */
public class RotamerMutator implements Mutator
{
    /** Not instantiable. */
    private RotamerMutator()
    {
        throw new IllegalArgumentException("not instantiable!");
    }

    /**
     * Sets the chi angles of a peptide residue to the specified values.
     * @param peptide the starting structure
     * @param residue the residue to be mutated
     * @param chis the chis (chi1, chi2, ..., chiN)
     * @return the rotated peptide
     */
    public static Peptide setChis(Peptide peptide, Residue residue, List<Double> chis)
    {
        // check the residue belongs in this peptide
        AminoAcid currentAminoAcid = residue.aminoAcid;
        Residue currentResidue = residue;
        List<Residue> sequence = peptide.sequence;
        int residueIndex = sequence.indexOf(currentResidue);
        if ( residueIndex == -1 )
            throw new IllegalArgumentException("residue does not belong to this peptide");
        return setChis(peptide, residueIndex, chis); 
    }

    /**
     * Sets the chi angles of a peptide residue to the specified values.
     * @param peptide the starting structure
     * @param residueIndex the index of the residue to adjust
     * @param chis the chis (chi1, chi2, ..., chiN)
     * @return the rotated peptide
     */
    public static Peptide setChis(Peptide peptide, int residueIndex, List<Double> chis)
    {
        // nulls are not allowed
        if ( peptide == null || chis == null )
            throw new NullPointerException("arguments cannot be null");
 
        // get necessary fields
        Residue residue = peptide.sequence.get(residueIndex);
        AminoAcid currentAminoAcid = residue.aminoAcid;
       
        // check the number of residues
        if ( chis.size() != residue.chis.size() )
            throw new IllegalArgumentException(String.format("chi list size mismatch (expected %d but found %d, %s)",
                                               residue.chis.size(), chis.size(), residue.description));

        // if the chis list is empty, merely return the original peptide
        if ( chis.size() == 0 )
            return peptide;

        Peptide newPeptide = peptide;
        LinkedList<ProtoTorsion> oldProtoTorsions = new LinkedList<>(residue.chis);
        LinkedList<IndexTorsion> indexTorsions = new LinkedList<>();
        LinkedList<Double> newChis = new LinkedList<>(chis);

        // make the mutation
        if ( currentAminoAcid.isProline() )
            {
                // disconnect bond in temporary molecule
                Molecule tempMolecule = peptide.moveAtoms(new HashMap<Atom,Atom>());
                ProtoTorsion chi3 = residue.chis.get(2);
                Atom atom3 = chi3.atom3;
                Atom atom4 = chi3.atom4;
                DefaultWeightedEdge e = tempMolecule.connectivity.removeEdge(atom3,atom4);
                if ( e == null )
                    throw new NullPointerException("unexpected null for proline connection");

                // get IndexTorsions
                indexTorsions.add(IndexTorsion.createIndexTorsion(oldProtoTorsions.get(0), tempMolecule));
                indexTorsions.add(IndexTorsion.createIndexTorsion(oldProtoTorsions.get(1), tempMolecule));

                // because proline contains a ring, special procedure
                // we don't set chi3 because it's hard and barely matters--would involve changing
                // bond lengths/angles as well and it just changes the pyrrolidine ring pucker
                oldProtoTorsions.removeLast();
                newChis.removeLast();

                // adjust the chis
                for (int j=0; j < indexTorsions.size(); j++)
                    {
                        IndexTorsion i = indexTorsions.get(j);
                        double newChiValue = newChis.get(j);
                        tempMolecule = tempMolecule.setDihedral(i, newChiValue);
                    }

                // move the atoms
                // note that this way of doing things preserves the previous connectivity
                newPeptide = peptide.setMolecule(tempMolecule);
            }
        else
            {
                // get IndexTorsions
                for (ProtoTorsion t : oldProtoTorsions)
                    indexTorsions.add(IndexTorsion.createIndexTorsion(t, peptide));

                // can proceed directly to dihedral angle setting
                Molecule tempMolecule = (Molecule)peptide;
                for (int j=0; j < newChis.size(); j++)
                    {
                        IndexTorsion i = indexTorsions.get(j);
                        double newTorsionValue = newChis.get(j);
                        tempMolecule = tempMolecule.setDihedral(i,newTorsionValue);
                    }

                newPeptide = peptide.setMolecule(tempMolecule);
            }

        // for debugging: print out old dihedral and new dihedrals
        //for (ProtoTorsion p : residue.chis)
        //    System.out.println("Old: " + p.getDihedralAngle());
        //for (ProtoTorsion p : newPeptide.sequence.get(residueIndex).chis)
        //    System.out.println("New: " + p.getDihedralAngle());

        // return the new peptide
        return newPeptide;
    }

    /**
     * Takes a residue and randomly alters its sidechain chi angles using the
     * Dunbrack rotamer library.
     * @param peptide the peptide to be mutated
     * @param i the index of the residue to be mutated
     * @return the mutated peptide
     */
    public static Peptide mutateChis(Peptide peptide, int i)
    {
        return mutateChis(peptide, peptide.sequence.get(i));
    }

    /**
     * Takes a residue and randomly alters its sidechain chi angles using the
     * Dunbrack rotamer library.
     * @param peptide the peptide to be mutated
     * @param residue the residue in peptide to be changed
     * @return the mutated peptide
     */
    public static Peptide mutateChis(Peptide peptide, Residue residue)
    {
        // get list of chis
        List<Double> newChis = null;
        double omegaValue = residue.omega.getDihedralAngle();
        double phiValue = residue.phi.getDihedralAngle();
        double psiValue = residue.psi.getDihedralAngle();

        if ( residue.aminoAcid.rotamerType == AminoAcid.RotamerType.HAS_NO_ROTAMERS )
            {
                // no rotamers, so don't do anything
                return peptide;
            }
        else if ( residue.aminoAcid == AminoAcid.TS )
            {
                // for transition states, just take a wild guess
                List<Double> specialChis = RotamerDatabase.getRandomRotamer(AminoAcid.SER, omegaValue, phiValue, psiValue);
                newChis = ImmutableList.of(specialChis.get(0), 120.0);
            }
        else if ( residue.aminoAcid == AminoAcid.DPRO ) 
            {
                // for D-proline
                // these numbers have already been inverted so un-invert them
                newChis = RotamerDatabase.getRandomRotamer(AminoAcid.LPRO, omegaValue, -1.0*phiValue, -1.0*psiValue);
                // invert the answer we get
                newChis = ImmutableList.of(newChis.get(0)*-1.0, newChis.get(1)*-1.0, newChis.get(2)*-1.0);
            }
        else
            {
                // for normal amino acids
                newChis = RotamerDatabase.getRandomRotamer(residue.aminoAcid, omegaValue, phiValue, psiValue);
            }

        if ( newChis.size() != residue.chis.size() )
            throw new IllegalArgumentException("chi list size mismatch! residue: " + residue.aminoAcid.shortName + " newChis: " + newChis.size() + " residue chis: " + residue.chis.size());
        
        // for debugging
        //for (Double d : newChis)
        //    System.out.println("requested: " + d);

        return setChis(peptide, residue, newChis);
    }
}
