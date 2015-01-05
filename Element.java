/**
 * This class represents an element.
 */
public enum Element
{
    CARBON  ("C"),
    HYDROGEN("H"),
    NITROGEN("N"),
    OXYGEN  ("O"),
    SULFUR  ("S");

    /** the atomic symbol for this element */
    public String symbol;

    /** Constructor. */
    Element(String symbol)
    {
        this.symbol = symbol;
    }

    /**
     * Returns the Element corresponding to the inputted String.
     * @param symbol the symbol for the requested element
     * @return the corresponding Element enum element
     */
    public static Element getElement(String symbol)
    {
        if ( symbol.equals("C") || symbol.equals("CA") )
            return CARBON;
        else if ( symbol.equals("H") || symbol.equals("HN") || symbol.equals("HS") || symbol.equals("HO") )
            return HYDROGEN;
        else if ( symbol.equals("N") )
            return NITROGEN;
        else if ( symbol.equals("O") || symbol.equals("O-") || symbol.equals("OH") )
            return OXYGEN;
        else if ( symbol.equals("S") || symbol.equals("SH") || symbol.equals("SS") )
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
