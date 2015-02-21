import java.util.*;
import java.io.*;
import com.google.common.collect.*;

public class EnergyBreakdown implements Immutable, Serializable, Result
{
    /** For serialization. */
    public static final long serialVersionUID = 1L;

    /** The total energies in kcal by residue, from the N- to C-terminus. */
    public final List<Double> energyByResidue;

    /** The total energy in kcal. */
    public final double totalEnergy;
    
    /** The total solvation energy in kcal. */
    public final double solvationEnergy;

    /** The total potential energy in kcal. */
    public final double potentialEnergy;

    /** The exposed surface area of each atom in the same order as in molecule.contents. */
    public final List<Double> sasa;

    /** The kind of energies represented here (OPLS or AMOEBA). */
    public final Forcefield type;

    /** Represents a blank energy breakdown. */
    public static final EnergyBreakdown BLANK = new EnergyBreakdown(null, 0.0, 0.0, 0.0, null, null); 

    public EnergyBreakdown(List<Double> energyByResidue, double totalEnergy, double solvationEnergy, double potentialEnergy, List<Double> sasa, Forcefield type)
    {
	    this.energyByResidue = energyByResidue;
        this.totalEnergy = totalEnergy;
        this.solvationEnergy = solvationEnergy;
        this.potentialEnergy = potentialEnergy;
        if ( totalEnergy != potentialEnergy + solvationEnergy )
            throw new IllegalArgumentException("energies don't add up");
        this.sasa = sasa;
        this.type = type;
    }

    /**
     * Creates a report of the reference energies for a peptide.
     * @param peptide the peptide this EnergyBreakdown was created from
     * @return the line by line listing of what all the reference energies are
     */
    public String createReportString(Peptide peptide)
    {
        if ( peptide.sequence.size() != energyByResidue.size() || this == BLANK )
            throw new IllegalArgumentException("wrong peptide!");

        String reportString = new Date().toString() + "\n\nReference Energies:\n";
        List<Residue> sequence = peptide.sequence;
        for (int i=0; i < energyByResidue.size(); i++)
            reportString = reportString + String.format("   %-30s %15.4f\n", sequence.get(i).description, energyByResidue.get(i));
        reportString = reportString + String.format("   Total Energy: %.4f\n", totalEnergy);
        return reportString;
    }

    public EnergyBreakdown addReferenceEnergy(double referenceEnergy)
    {
        return new EnergyBreakdown(null, totalEnergy-referenceEnergy, 0.0, totalEnergy-referenceEnergy, null, type);
    }

    @Override
    public int hashCode()
    {
        return Objects.hash(energyByResidue, totalEnergy, solvationEnergy, potentialEnergy, sasa, type);
    }

    public boolean equals(Object obj)
    {
        if ( obj == null )
            return false;
        if ( obj == this )
            return true;
        if ( !(obj instanceof EnergyBreakdown) )
            return false;

        EnergyBreakdown b = (EnergyBreakdown)obj;
        if ( Objects.equals(this.energyByResidue, b.energyByResidue) &&
                            this.totalEnergy == b.totalEnergy &&
                            this.solvationEnergy == b.solvationEnergy &&
                            this.potentialEnergy == b.potentialEnergy &&
             Objects.equals(this.sasa, b.sasa) &&
                            this.type == b.type                          )
            return true;
        return false;
    }

    @Override
    public String toString()
    {
	    return energyByResidue.toString() + "\nTotal Energy: " + totalEnergy;
    }
}
