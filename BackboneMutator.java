import org.apache.commons.math3.geometry.euclidean.threed.*;
import java.io.*;
import java.util.*;
import org.apache.commons.math3.distribution.NormalDistribution;
import java.util.concurrent.*;
import org.jgrapht.*;
import org.jgrapht.graph.*;
import org.jgrapht.alg.*;
import com.google.common.collect.*;

/**
 * This collects together methods for taking a peptide and changing its
 * omega, phi, and psi angles to appropriate randomly-selected values.
 */
public class BackboneMutator implements Mutator
{
    /** This is the probability for giving a cis omega dihedral angle for a proline. */
    public static final double CIS_PROLINE_PROBABILITY = 0.08;

    /** This is the width of the normal distribution in degrees to draw omega angles for if no database data are available. */
    public static final double DEFAULT_OMEGA_WIDTH = 5.0;

    /** should not be invoked */
    private BackboneMutator()
    {
        throw new IllegalArgumentException("not instantiable");
    }

    /**
     * Takes a peptide and sets the omega angle.  We don't check that residue belongs to peptide.
     * @param peptide the peptide to set the omega of
     * @param residue the residue to set the omega of
     * @param omegaAngle the angle to set the omega to in degrees
     * @return the rotated peptide
     */
    public static Peptide setOmega(Peptide peptide, Residue residue, double omegaAngle)
    {
        ProtoTorsion omega = residue.omega;
        return peptide.setMolecule( peptide.setDihedral(omega, omegaAngle) );
    }

    /**
     * Takes a peptide and sets the omega angle.  We don't check that residue belongs to peptide.
     * @param peptide the peptide to set the omega of
     * @param i the index of the residue to adjust
     * @param omegaAngle the angle to set the omega to in degrees
     * @return the rotated peptide
     */
    public static Peptide setOmega(Peptide peptide, int i, double omegaAngle)
    {
        ProtoTorsion omega = peptide.sequence.get(i).omega;
        return peptide.setMolecule( peptide.setDihedral(omega, omegaAngle) );
    }

    /**
     * Convenience method that mutates the omega torsion angle of the i-th
     * residue (i=0,1,...,n-1).  The omega angle is randomly selected.
     * @param peptide the peptide to mutate
     * @param i the index of the residue to mutate
     * @return the mutated peptide
     */
    public static Peptide mutateOmega(Peptide peptide, int i)
    {
        if ( i < 0 || i > peptide.sequence.size() - 1 )
            throw new IllegalArgumentException("residue index out of range");
        AminoAcid aa2 = peptide.sequence.get(i).aminoAcid;

        // get the residue to the left
        Double omegaValue = null;
        if ( i > 0 )
            {
                if ( aa2.isProline() )
                    {
                        // give this some chance of becoming cis-proline
                        double draw = ThreadLocalRandom.current().nextDouble();
                        if ( draw < CIS_PROLINE_PROBABILITY )
                            omegaValue = new NormalDistribution(0.0, DEFAULT_OMEGA_WIDTH).sample();
                    }
                if ( omegaValue == null )
                    {
                        // this is a normal amino acid, so use information from OmegaDatabase
                        AminoAcid aa1 = peptide.sequence.get(i-1).aminoAcid;
                        double psi0   = peptide.sequence.get(i-1).psi.getDihedralAngle();
                        double phi1   = peptide.sequence.get(i).phi.getDihedralAngle();
                        omegaValue    = OmegaDatabase.getOmega(aa1, psi0, aa2, phi1);
                    }
            }
        else
            {
                // there is no residue to the left, so just return a normally distributed value with a small width
                omegaValue = new NormalDistribution(180.0, DEFAULT_OMEGA_WIDTH).sample();
            }
        
        // make the change
        return setOmega(peptide, peptide.sequence.get(i), omegaValue);
    }

    /**
     * Takes a peptide and sets the (phi,psi) value to a likely value using Dunbrack Ramachandran data.
     * @param peptide the starting structure
     * @param i the index of the residue to be mutated
     * @return the result of the mutation
     */
    public static Peptide mutatePhiPsi(Peptide peptide, int i)
    {
        return mutatePhiPsi(peptide, peptide.sequence.get(i));
    }

    /**
     * Takes a peptide and sets the (phi,psi) value to a likely value using Dunbrack Ramachandran data.
     * @param peptide the starting structure
     * @param residue the residue whose (phi,psi) angles are to be mutated
     * @return the result of the mutation
     */
    public static Peptide mutatePhiPsi(Peptide peptide, Residue residue)
    {
        // get residue index
        List<Residue> sequence = peptide.sequence;
        int residueIndex = sequence.indexOf(residue);

        // get adjacent amino acids
        AminoAcid centralAminoAcid = residue.aminoAcid;
        if ( centralAminoAcid == AminoAcid.DPRO )
            centralAminoAcid = AminoAcid.LPRO;
        double centralOmega = residue.omega.getDihedralAngle();

        Double newPhiValue = null;
        Double newPsiValue = null;
        
        AminoAcid leftAminoAcid = null;
        Double leftOmega = null;
        if ( residueIndex > 0 )
            {
                Residue leftResidue = sequence.get(residueIndex-1);
                leftAminoAcid = leftResidue.aminoAcid;
                leftOmega = leftResidue.omega.getDihedralAngle();
            }

        AminoAcid rightAminoAcid = null;
        Double rightOmega = null;
        if ( residueIndex < sequence.size() - 1 )
            {
                Residue rightResidue = sequence.get(residueIndex+1);
                rightAminoAcid = rightResidue.aminoAcid;
                rightOmega = rightResidue.omega.getDihedralAngle();
            }

        // if D-proline is adjacent, set to null because we don't have data for that
        if ( leftAminoAcid == AminoAcid.DPRO )
            leftAminoAcid = null;
        if ( rightAminoAcid == AminoAcid.DPRO )
            rightAminoAcid = null;

        // use serine data for adjacent transition states
        if ( leftAminoAcid == AminoAcid.TS )
            leftAminoAcid = AminoAcid.SER;
        if ( rightAminoAcid == AminoAcid.TS )
            rightAminoAcid = AminoAcid.SER;

        // obtain appropriate probability distribution of phi,psi
        DiscreteProbabilityDistribution<RotamerLibrary.Angles> DPD = null;
        if ( leftAminoAcid != null && rightAminoAcid != null )
            DPD = RamachandranDatabase.getTripletDistribution(leftAminoAcid, leftOmega, centralAminoAcid, centralOmega, rightAminoAcid, rightOmega);
        else if ( leftAminoAcid == null && rightAminoAcid != null )
            DPD = RamachandranDatabase.getRightDistribution(centralAminoAcid, centralOmega, rightAminoAcid, rightOmega);
        else if ( leftAminoAcid != null && rightAminoAcid == null )
            DPD = RamachandranDatabase.getLeftDistribution(leftAminoAcid, leftOmega, centralAminoAcid, centralOmega);
        else
            throw new IllegalArgumentException("should not be possible");
    
        // draw random phi,psi from the distribution
        RotamerLibrary.Angles draw = DPD.getRandom();
        newPhiValue = draw.phi;
        newPsiValue = draw.psi;
        
        if ( residue.aminoAcid == AminoAcid.DPRO )
            {
                // no rotamer data for non-hairpin D-proline so get data for L-proline in this position and invert numbers
                newPhiValue = newPhiValue * -1.0;
                newPsiValue = newPsiValue * -1.0;
                //System.out.println(String.format("phi,psi for dpro: %6.1f, %6.1f", newPhiValue, newPsiValue));
            }

        // make the mutation
        return setPhiPsi(peptide, residue, newPhiValue, newPsiValue);
    }

    /**
     * Sets the backbone (phi,psi) angles of a residue to specific values.
     * @param peptide the starting structure
     * @param residue the residue whose (phi,psi) angles are to be mutated
     * @param newPhiValue the phi value we will be setting residue's phi to in degrees
     * @param newPsiValue the psi value we will be setting residue's psi to in degeres
     * @return the result of the mutation
     */
    public static Peptide setPhiPsi(Peptide peptide, Residue residue, double newPhiValue, double newPsiValue)
    {
        int residueIndex = peptide.sequence.indexOf(residue);
        if ( residueIndex == -1 )
            throw new IllegalArgumentException("this residue does not belong to the specified peptide");
        
        return setPhiPsi(peptide, residueIndex, newPhiValue, newPsiValue);
    }

    /**
     * Sets the backbone (phi,psi) angles of a residue to specific values.
     * @param peptide the starting structure
     * @param residueIndex the index of the residue whose (phi,psi) angles are to be mutated
     * @param newPhiValue the phi value we will be setting residue's phi to in degrees
     * @param newPsiValue the psi value we will be setting residue's psi to in degeres
     * @return the result of the mutation
     */
    public static Peptide setPhiPsi(Peptide peptide, int residueIndex, double newPhiValue, double newPsiValue)
    {
        if ( newPhiValue < -180.0 || newPhiValue > 180.0 || newPsiValue < -180.0 || newPsiValue > 180.0 )
            throw new IllegalArgumentException("angle out of range");
        
        Residue residue = peptide.sequence.get(residueIndex);
        ProtoTorsion phi = residue.phi;
        ProtoTorsion psi = residue.psi;

        double omegaValue = residue.omega.getDihedralAngle();

        Peptide newPeptide = null;
        if ( residue.aminoAcid.isProline() )
            {
                // special procedure for proline
                //System.out.println("old chi1: " + residue.chis.get(0).getDihedralAngle());
                //System.out.println("old chi2: " + residue.chis.get(1).getDihedralAngle());
                //System.out.println("old chi3: " + residue.chis.get(2).getDihedralAngle());

                // disconnect bond in temporary molecule
                Molecule tempMolecule = peptide.moveAtoms(new HashMap<Atom,Atom>());
                ProtoTorsion chi3 = residue.chis.get(2);
                Atom atom3 = chi3.atom3;
                Atom atom4 = chi3.atom4;
                DefaultWeightedEdge e = tempMolecule.connectivity.removeEdge(atom3,atom4);
                if ( e == null )
                    throw new NullPointerException("unexpected null for proline connection");

                // keep psi for later
                IndexTorsion psiIndexTorsion = IndexTorsion.createIndexTorsion(psi, tempMolecule);
                
                // save information for chi mutations
                LinkedList<ProtoTorsion> oldProtoTorsions = new LinkedList<>(residue.chis);
                oldProtoTorsions.removeLast();
                List<IndexTorsion> indexTorsions = new LinkedList<>();
                for (ProtoTorsion p : oldProtoTorsions)
                    indexTorsions.add(IndexTorsion.createIndexTorsion(p, tempMolecule));

                // make phi mutation
                tempMolecule = tempMolecule.setDihedral(phi, newPhiValue);

                // make psi mutation
                tempMolecule = tempMolecule.setDihedral(psiIndexTorsion, newPsiValue);

                // choose new proline chis
                // there's no easy way to adjust the third chi value, so we don't set it
                List<Double> newChis = null;
                if ( residue.aminoAcid == AminoAcid.LPRO )
                    {
                        newChis = RotamerDatabase.getRandomRotamer(AminoAcid.LPRO, omegaValue, newPhiValue, newPsiValue);;
                        newChis = ImmutableList.of(newChis.get(0), newChis.get(1));
                    }
                else if ( residue.aminoAcid == AminoAcid.DPRO )
                    {
                        // these numbers have already been inverted so un-invert them
                        newChis = RotamerDatabase.getRandomRotamer(AminoAcid.LPRO, omegaValue, -1.0*newPhiValue, -1.0*newPsiValue);;
                        // invert the answer we get
                        newChis = ImmutableList.of(newChis.get(0)*-1.0, newChis.get(1)*-1.0);
                    }
                else
                    throw new IllegalArgumentException("I don't know how to deal with this kind of proline");

                // adjust the chis
                for (int j=0; j < oldProtoTorsions.size(); j++)
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
                // make mutation directly
                IndexTorsion psiIndexTorsion = IndexTorsion.createIndexTorsion(psi, peptide);
                newPeptide = peptide.setMolecule( peptide.setDihedral(phi, newPhiValue) );
                Residue newResidue = newPeptide.sequence.get(residueIndex);
                newPeptide = newPeptide.setMolecule( newPeptide.setDihedral(psiIndexTorsion, newPsiValue) );
            }
        return newPeptide;
    }
}
