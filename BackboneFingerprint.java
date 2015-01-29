import java.util.*;
import com.google.common.collect.*;

/**
 * This immutable class can be used as a key to remove duplicates from a collection of peptide backbones.
 */
public class BackboneFingerprint implements Immutable
{
    /** How much to round (phi,psi) backbone angles to for duplicate elimination.  In degrees. */
    public static final double PHI_PSI_ANGLE_TOLERANCE = 10.0;

    /** How much to round (phi,psi) backbone angles to for duplicate elimination.  In degrees. */
    public static final double OMEGA_ANGLE_TOLERANCE = 5.0;

    /** The omega angle for every residue. */
    public final List<Integer> omegas;
    
    /** The phi angle for every residue. */
    public final List<Integer> phis;

    /** The psi angle for every residue. */
    public final List<Integer> psis;

    /** 
     * Create a BackboneFingerprint from a peptide.
     * Used to create an initial BackboneFingerprint that will be changed by HaripinJobs
     * @param peptide the input peptide for which to generate a backbone fingerprint
     */
    public BackboneFingerprint(Peptide peptide)
    {
        List<Integer> tempOmegas = new ArrayList<>();
        List<Integer> tempPhis = new ArrayList<>();
        List<Integer> tempPsis = new ArrayList<>();

        for (Residue r : peptide.sequence)
        {
            int omega = round(r.omega.getDihedralAngle(), OMEGA_ANGLE_TOLERANCE);
            int phi = round(r.phi.getDihedralAngle(), PHI_PSI_ANGLE_TOLERANCE);
            int psi = round(r.psi.getDihedralAngle(), PHI_PSI_ANGLE_TOLERANCE);
            tempOmegas.add(omega);
            tempPhis.add(phi);
            tempPsis.add(psi);
        }

        for (Integer d : tempOmegas)
            if ( d < -180.0 || d > 180.0 )
                throw new IllegalArgumentException("unexpected omega: " + tempOmegas.toString());
        for (Integer d : tempPsis)
            if ( d < -180.0 || d > 180.0 )
                throw new IllegalArgumentException("unexpected psi: " + tempPsis.toString());
        for (Integer d : tempPhis)
            if ( d < -180.0 || d > 180.0 )
                throw new IllegalArgumentException("unexpected omega: " + tempPhis.toString());
               
        this.omegas = ImmutableList.copyOf(tempOmegas);
        this.phis = ImmutableList.copyOf(tempPhis);
        this.psis = ImmutableList.copyOf(tempPsis);
    }

    /** 
     * Modify a set of backbone angles in BackboneFingerprint.
     * Will return a new BackboneFingerprint because this class is Immutable
     * @param residueIndex the index of the residue whose backbone angles will be changed
     * @return a BackboneFingerprint with updated backbone angles from the passed in values
     */
    public BackboneFingerprint changeResidueBackboneAngles(int residueIndex, double omega, double phi, double psi)
    {
        if (omega < -180.0 || omega > 180.0)
            throw new IllegalArgumentException("Invalid omega");
        if (phi < -180.0 || phi > 180.0)
            throw new IllegalArgumentException("Invalid phi");
        if (psi < -180.0 || psi > 180.0)
            throw new IllegalArgumentException("Invalid psi");
           
        List<Integer> tempOmegas = new ArrayList<>(omegas);
        List<Integer> tempPhis = new ArrayList<>(phis);
        List<Integer> tempPsis = new ArrayList<>(psis);

        tempOmegas.set(residueIndex, round(omega, OMEGA_ANGLE_TOLERANCE));
        tempPhis.set(residueIndex, round(phi, PHI_PSI_ANGLE_TOLERANCE));
        tempPsis.set(residueIndex, round(psi, PHI_PSI_ANGLE_TOLERANCE));

        return new BackboneFingerprint(tempOmegas, tempPhis, tempPsis);
    }

    /**
     * Create a BackboneFingerprint by calculating the backbone angles of the active residues.
     */
    public BackboneFingerprint(List<Integer> omegas, List<Integer> phis, List<Integer> psis)
    {
        if ( omegas == null || omegas.size() == 0 )
            throw new NullPointerException("empty angles"); 
        if ( phis == null || phis.size() == 0 )
            throw new NullPointerException("empty angles"); 
        if ( psis == null || psis.size() == 0 )
            throw new NullPointerException("empty angles"); 
                    
        this.omegas = ImmutableList.copyOf(omegas);
        this.phis = ImmutableList.copyOf(phis);
        this.psis = ImmutableList.copyOf(psis);
    }

    /**
     * Round a number to the nearest tolerance.
     * @param d the number to round
     * @param tolerance round the number to this tolerance (degrees)
     * @return the rounded number
     */
    private static int round(double d, double tolerance)
    {
        return (int) (Math.round(d/tolerance) * tolerance);
    }

    @Override
    public String toString()
    {
        return "Omegas: " + omegas.toString() + "\nPhis: " + phis.toString() + "\nPsis: " + psis.toString();
    }

    @Override
    public int hashCode()
    {
        return Objects.hash(phis, psis, omegas);
    }

    @Override
    public boolean equals(Object obj)
    {
        if ( obj == null )
            return false;
        if ( obj == this ) 
            return true;
        if ( !(obj instanceof BackboneFingerprint) )
            return false;

        BackboneFingerprint f = (BackboneFingerprint)obj;
        if ( omegas.equals(f.omegas) && phis.equals(f.phis) && psis.equals(f.psis))
            return true;
        return false;
    }
}
