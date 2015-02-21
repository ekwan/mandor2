import java.util.*;
import com.google.common.collect.*;

/**
 * This enum represents an element.
 */
public enum Element
{
    CARBON  ("C"),
    HYDROGEN("H"),
    NITROGEN("N"),
    OXYGEN  ("O"),
    SULFUR  ("S");

    public static final Set<String> CARBON_STRINGS   = ImmutableSet.of("C", "CA", "CB", "CT", "CM", "CN", "CO", "C*", "CP", "CV", "CW");
    public static final Set<String> HYDROGEN_STRINGS = ImmutableSet.of("H", "HA", "HN", "HS", "HO", "HC", "H2", "H3");
    public static final Set<String> NITROGEN_STRINGS = ImmutableSet.of("N", "N2", "N3", "NB", "NA");
    public static final Set<String> OXYGEN_STRINGS   = ImmutableSet.of("O", "O-", "OH", "OS", "O2");
    public static final Set<String> SULFUR_STRINGS   = ImmutableSet.of("S", "SH", "SS");

    /** The atomic symbol. */
    public final String symbol;

    /**
     * Constructs an element.
     * @param symbol the atomic symbol
     */
    Element(String symbol)
    {
        this.symbol = symbol;
    }

    /**
     * Identifies the element corresponding to a string.
     * @param symbol the symbol for the requested element
     * @return the corresponding Element enum element
     */
    public static Element getElement(String symbol)
    {
        if ( CARBON_STRINGS.contains(symbol) )
            return CARBON;
        else if ( HYDROGEN_STRINGS.contains(symbol) )
            return HYDROGEN;
        else if ( NITROGEN_STRINGS.contains(symbol) )
            return NITROGEN;
        else if ( OXYGEN_STRINGS.contains(symbol) )
            return OXYGEN;
        else if ( SULFUR_STRINGS.contains(symbol) )
            return SULFUR;
        else
            throw new IllegalArgumentException("Unrecognized atom symbol (" + symbol + ")!");
    }

    /**
     * Returns a string representation of this Element.
     * @return the String
     */
    public String toString()
    {
        return String.format("%s", symbol);
    }
}
