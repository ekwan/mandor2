import java.util.*;
import java.util.concurrent.*;
import com.google.common.collect.*;
import java.io.*;
import com.google.common.collect.*;
import org.apache.commons.math3.geometry.euclidean.threed.*;

/** 
 * This class is meant as a key to eliminate duplicates in collections of Peptide conformations.
 * All omega, phi, psi, and chi angles are stored and used for geometry comparison.  The energy
 * can also be stored but it is not used in the comparison.
 */
public class PeptideFingerprint implements Serializable
{
    /** For serialization. */
    public static final long serialVersionUID = 1L;

    /** Maximum amount in degrees two dihedral angles can be different for two PoseEnergyResults to be considered the same. */
    public static final double ANGLE_TOLERANCE = 15.0;

    /** For comparing the geometries of two peptides. */
    private final short[] geometrySignature;

    /** The energy of the peptide. */
    public final double energy;

    /**
     * Generates a Pose Fingerprint from a template peptide, a List<Vector3D> of the pose's atom position, and an energy
     * @param templatePeptide the peptide that is the template for structure
     * @param structure a conformer of the template peptide
     * @param energy the pose energy which is defined as the sum of the Tinker and Omnisol solvation energies
     */
    public PeptideFingerprint(Peptide templatePeptide, List<Vector3D> structure, double energy)
    {
        this.energy = energy;
        if ( templatePeptide.contents.size() != structure.size() )
            throw new IllegalArgumentException("size mismatch");
        this.geometrySignature = getGeometrySignature(templatePeptide, structure);
    }
    
    /**
     * Generates a Pose Fingerprint from the a peptide object with the pose and an energy
     * @param pose a pose of the ground or transition state
     * @param energy the pose energy which is defined as the sum of the Tinker and Omnisol solvation energies
     */
    public PeptideFingerprint(Peptide pose, double energy)
    {
        this.energy = energy;
        this.geometrySignature = getGeometrySignature(pose);
    }
    
    /**
     * Rounds the given angle to the specified tolerance.
     * @param torsion the torsion to get the dihedral angle for
     * @return the torsion angle rounded to the nearest {@link #ANGLE_TOLERANCE}
     */
    public static short getRoundedAngle(ProtoTorsion torsion)
    {
        double angle = torsion.getDihedralAngle();
        double rounded = Math.rint(angle/ANGLE_TOLERANCE)*ANGLE_TOLERANCE + 0.0;
        return (short)rounded;
    }

    /**
     * Rounds the given angle to the specified tolerance.
     * @param torsion the torsion to get the dihedral angle for
     * @return the torsion angle rounded to the nearest {@link #ANGLE_TOLERANCE}
     */
    public static short getRoundedAngle(Peptide peptide, ProtoTorsion torsion, List<Vector3D> structure)
    {
        int index1 = peptide.contents.indexOf(torsion.atom1);
        int index2 = peptide.contents.indexOf(torsion.atom2);
        int index3 = peptide.contents.indexOf(torsion.atom3);
        int index4 = peptide.contents.indexOf(torsion.atom4);
        double angle = AbstractTorsion.getDihedralAngle(structure.get(index1), structure.get(index2), structure.get(index3), structure.get(index4));
        double rounded = Math.rint(angle/ANGLE_TOLERANCE)*ANGLE_TOLERANCE + 0.0;
        return (short)rounded;
    }

    /**
     * Creates a list of doubles that contains the backbone and sidechain dihedral angles.
     * @param peptide the peptide to analyze
     * @return the list of angles: blocks of phi, psi, chi1, chi2, ... from residue 0 onwards
     */
    public static short[] getGeometrySignature(Peptide peptide)
    {
        List<Short> angles = new ArrayList<>(peptide.sequence.size()*5);
        for (Residue r : peptide.sequence)
            {
                angles.add(getRoundedAngle(r.omega));
                angles.add(getRoundedAngle(r.phi));
                angles.add(getRoundedAngle(r.psi));
                for (ProtoTorsion chi : r.chis)
                    angles.add(getRoundedAngle(chi));
            }
        short[] returnArray = new short[angles.size()];
        for (int i=0; i < angles.size(); i++)
            returnArray[i] = (short)angles.get(i);
        return returnArray;
    }

    /**
     * Creates a list of doubles that contains the backbone and sidechain dihedral angles.
     * @param templatePeptide the peptide with its contents parallel to structure
     * @param structure a pose typically from Ground State
     * @return the list of angles: blocks of phi, psi, chi1, chi2, ... from residue 0 onwards
     */
    private static short[] getGeometrySignature(Peptide templatePeptide, List<Vector3D> structure)
    {
        List<Short> angles = new ArrayList<>(templatePeptide.sequence.size() * 5);
        for (Residue r : templatePeptide.sequence)
        {
            angles.add(getRoundedAngle(templatePeptide, r.omega, structure));
            angles.add(getRoundedAngle(templatePeptide, r.phi, structure));
            angles.add(getRoundedAngle(templatePeptide, r.psi, structure));
            for (ProtoTorsion chi : r.chis)
                angles.add(getRoundedAngle(templatePeptide, chi, structure));
        }
        short[] returnArray = new short[angles.size()];
        for (int i=0; i < angles.size(); i++)
            returnArray[i] = (short)angles.get(i);
        return returnArray;
    }
    
    @Override 
    public boolean equals(Object obj)
    {
        if ( obj == null)
            return false;
        if ( obj == this)
            return false;
        if ( !(obj instanceof PeptideFingerprint))
            return false;
        PeptideFingerprint r = (PeptideFingerprint)obj;

        if ( geometrySignature.equals(r.geometrySignature) )
            return true;
        return false;
    }

    /**
     * Returns the hash code.
     * @return the hash code
     */
    @Override
    public int hashCode()
    {
        return Objects.hash(geometrySignature);
    }


    @Override
    public String toString()
    {
        return String.format("%s : %.2f", Arrays.toString(geometrySignature), energy);
    }
}
