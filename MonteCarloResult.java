import java.util.*;
import com.google.common.collect.*;

public class MonteCarloResult implements Result, Immutable
{
    /** the best results from the minimization */
    public final ImmutableList<Peptide> bestPeptides;
    
    public MonteCarloResult(List<Peptide> bestPeptides)
    {
        // check that each peptide has an EnergyBreakdown
        for (Peptide p : bestPeptides)
            if ( p.energyBreakdown == null || p.energyBreakdown == EnergyBreakdown.BLANK )
                throw new IllegalArgumentException("missing energy breakdown for peptide " + p.name);
        this.bestPeptides = ImmutableList.copyOf(bestPeptides);
    }

    @Override
    public int hashCode()
    {
        return Objects.hash(bestPeptides);
    }

    @Override
    public boolean equals(Object obj)
    {
        if ( obj == null )
            return false;
        if ( obj == this )
            return true;
        if ( !(obj instanceof MonteCarloResult) )
            return false;

        MonteCarloResult r = (MonteCarloResult)obj;
        if ( this.bestPeptides.equals(r.bestPeptides) )
            return true;
        return false;
    }

    @Override
    public String toString()
    {
        String returnString = "[";
        for (Peptide p : bestPeptides)
            returnString += String.format("%.1f, ", p.energyBreakdown.totalEnergy);
        if ( returnString.length() > 1 )
            returnString = returnString.substring(0,returnString.length()-2);
        returnString += "]";
        return returnString;
    }
}
