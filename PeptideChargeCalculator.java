import java.util.*;
import com.google.common.collect.*;

/**
 * This is a utility class that calculates the charge of a peptide.
 */
public class PeptideChargeCalculator implements Immutable
{
    /** Residues containing these strings in their descriptions are positively charged by one unit. */
    public static Set<String> POSITIVELY_CHARGED = ImmutableSet.of("arginine", "lysine");

    /** Residues containing these strings in their descriptions are negatively charged by one unit. */
    public static Set<String> NEGATIVELY_CHARGED = ImmutableSet.of("aspartate", "glutamate");

    /**
     * Calculates the formal charge by looking at the names of the residues.
     * @param descriptions the descriptions to analyze
     * @return the overall charge
     */
    public static int getCharge(List<String> descriptions)
    {
        int totalCharge = 0;
        for (String description : descriptions)
            {
                for (String s : POSITIVELY_CHARGED)
                    {
                        if ( description.indexOf(s) > -1 )
                            {
                                totalCharge += 1;
                                break;
                            }
                    }
                for (String s : NEGATIVELY_CHARGED)
                    {
                        if ( description.indexOf(s) > -1 ) 
                            {
                                totalCharge += -1;
                                break;
                            }
                    }
            }
        return totalCharge;
    }
    
    /**
     * Calculates the formal charge of a peptide that is intended for rotamer reconstitution.
     * @param peptide the template
     * @param rotamers the rotamers that will be used on the template
     * @return the formal charge
     */
    public static int getCharge(Peptide peptide, List<Rotamer> rotamers)
    {
        List<String> descriptions = new ArrayList<>(peptide.sequence.size());
        for (Residue r : peptide.sequence)
            descriptions.add(r.description);
        for (Rotamer r : rotamers)
            descriptions.set(r.sequenceIndex, r.description);
        return getCharge(descriptions);
    }

    /**
     * Calculates the formal charge of a peptide.
     * @param peptide the peptide to analyze
     * @return the formal charge
     */
    public static int getCharge(Peptide peptide)
    {
        return getCharge(peptide, -1);
    }

    /**
     * Calculates the formal charge of a peptide, excluding some position.
     */
    public static int getCharge(Peptide peptide, int ignoreIndex)
    {
        List<String> descriptions = new ArrayList<>(peptide.sequence.size());
        for (int i=0; i < peptide.sequence.size(); i++)
            {
                if ( i == ignoreIndex )
                    continue;
                Residue r = peptide.sequence.get(i);
                descriptions.add(r.description);
            }
        return getCharge(descriptions);
    }
}
