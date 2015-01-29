import java.util.*;
import java.io.*;
import com.google.common.collect.*;
import org.apache.commons.math3.geometry.euclidean.threed.*;

/**
 * This class represents something that can calculate the solvent-exposed surface area of a molecule.
 */
public abstract class SurfaceAreaCalculator
{
    /**
     * Calculates the solvent accessible surface area of the specified molecule.
     * @param molecule the molecule to analyze
     * @return the SASA by atom in angstroms^2
     */
    public abstract List<Double> calculateSASA(Molecule molecule);

    /**
     * Returns the atomic radii for the given molecule.
     * @param molecule the molecule to analyze
     * @return the radii in angstroms ordered by atom index
     */
    public abstract List<Double> getRadii(Molecule molecule);
}
