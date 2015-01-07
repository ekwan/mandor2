import java.util.*;

/**
 * Represents backbone-dependent rotamer data for an amino acid.
 */
public abstract class RotamerLibrary
{
    /** Represents (phi,psi) backbone angles. */
    public static class Angles
    {
        public final double phi, psi;

        public Angles(Double phi, Double psi)
        {
            if ( phi < -180.0 || phi > 180.0 || psi < -180.0 || psi > 180.0 )
                throw new IllegalArgumentException("angle out of range");
            this.phi = phi;
            this.psi = psi;
        }

        @Override
        public String toString()
        {
            return String.format("%.0f, %.0f", phi, psi);
        }

        @Override
        public int hashCode()
        {
            return Objects.hash(phi,psi);
        }
        
        @Override
        public boolean equals(Object obj)
        {
            if ( obj == null )
                return false;
            if ( obj == this )
                return true;
            if ( !(obj instanceof Angles) )
                return false;

            Angles a = (Angles)obj;
            if ( phi == a.phi && psi == a.psi )
                return true;
            return false;
        }
    }   
}
