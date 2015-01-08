import java.util.*;
import java.io.*;
import com.google.common.collect.*;
import java.util.concurrent.*;

/**
 * Represents the chemical definition of an amino acid.  
 * Cis- and trans-proline are represented together, because they are chemically the same.
 * However, L- and D-proline, have different representations because they are chemically distinct.
 *
 * Backbone-dependent rotamer data will be read elsewhere from the Dunbrack library.
 * We use traditional rotamers for amino acids containing all sp3-sp3 bonds in their sidechains.
 * That is, we use ordered tuples (X<sub>1</sub>, X<sub>2</sub>, ..., X<sub>n</sub>) to represent
 * the backbone torsion angles, where n is the number of sidechain torsions.  (OHs
 * are not counted in the number of torsions.  We use a combination of an ordered
 * tuple and a DiscreteProbabilityDistribution to represent the non-rotameric
 * amino acids.  These amino acids contain an sp<sup>3</sup>-sp<sup>2</sup>
 * torsion at the end (e.g., phenylalanine). That is, we use an ordered tuple
 * (X<sub>1</sub>, X<sub>2</sub>, ..., X<sub>n-1</sub>) to represent the standard
 * rotamer part of the sidechain and then a probability distribution to represent
 * the terminal torsion.  Some amino acids do not contain any rotable bonds (e.g.
 * glycine) and therefore do not have associated library data.
 */
public enum AminoAcid
{
    ALA("Ala",             "alanine",              RotamerType.HAS_NO_ROTAMERS,     Chirality.L),
    GLY("Gly",             "glycine",              RotamerType.HAS_NO_ROTAMERS,     Chirality.ACHIRAL),
    VAL("Val",             "valine",               RotamerType.IS_ROTAMERIC,        Chirality.L),
    LEU("Leu",             "leucine",              RotamerType.IS_ROTAMERIC,        Chirality.L),
    ILE("Ile",             "isoleucine",           RotamerType.IS_ROTAMERIC,        Chirality.L),

    LPRO("L-Pro",          "L-proline",            RotamerType.IS_ROTAMERIC,        Chirality.L),
    DPRO("D-Pro",          "D-proline",            RotamerType.SPECIAL,             Chirality.D),

    PHE("Phe",             "phenylalanine",        RotamerType.NON_ROTAMERIC,       Chirality.L), 
    TYR("Tyr",             "tyrosine",             RotamerType.NON_ROTAMERIC,       Chirality.L),
    TRP("Trp",             "tryptophan",           RotamerType.NON_ROTAMERIC,       Chirality.L),
    SER("Ser",             "serine",               RotamerType.IS_ROTAMERIC,        Chirality.L),
    THR("Thr",             "threonine",            RotamerType.IS_ROTAMERIC,        Chirality.L),
    CYS("Cys",             "cysteine",             RotamerType.IS_ROTAMERIC,        Chirality.L),
    MET("Met",             "methionine",           RotamerType.IS_ROTAMERIC,        Chirality.L),
    ASN("Asn",             "aspargine",            RotamerType.NON_ROTAMERIC,       Chirality.L),
    GLN("Gln",             "glutamine",            RotamerType.NON_ROTAMERIC,       Chirality.L),
    LYS("Lys",             "lysine",               RotamerType.IS_ROTAMERIC,        Chirality.L),
    ARG("Arg",             "arginine",             RotamerType.IS_ROTAMERIC,        Chirality.L),
    HIS("His",             "histidine",            RotamerType.NON_ROTAMERIC,       Chirality.L),
    ASP("Asp",             "aspartate",            RotamerType.NON_ROTAMERIC,       Chirality.L),
    GLU("Glu",             "glutamate",            RotamerType.NON_ROTAMERIC,       Chirality.L),
    
    TS ("TS",              "transition state",     RotamerType.SPECIAL,             Chirality.L);

    /** An abbreviation like "Ala". */
    public final String shortName;

    /** A full name like "alanine". */
    public final String fullName;

    /** Indicates the kinds of sidechain torsions present. */
    public final RotamerType rotamerType;

    /** The chirality, like L or D. */
    public final Chirality chirality;

    /** Standard constructor. */
    AminoAcid(String shortName, String fullName, RotamerType rotamerType, Chirality chirality)
    {
        this.fullName = fullName;
        this.shortName = shortName;
        this.rotamerType = rotamerType;
        this.chirality = chirality;
    }

    /**
     * Identifies the amino acid corresponding to a string.
	 * @param input that contains a string with an amino acid name
	 * @return the amino acid that corresponds to the string input
    */
    public static AminoAcid getAminoAcid(String input)
    {
        return valueOf(input.toUpperCase());
    }

    @Override
    public String toString()
    {
        return shortName;
    }

    /**
     * Indicates whether the amino acid can be represented by standard rotamers,
     * is non-rotameric, or contains no rotatable bonds at all.  For a full description,
     * see the <a href="http://dunbrack.fccc.edu/bbdep2010/">Dunbrack backbone-dependent rotamer library</a> page.
     */
    public enum RotamerType
    {
        /**
         * Represents an amino acid that has standard rotameric degrees of freedom.
         * That is, all sidechain torsions involves sp3-sp3 bonds.
         */
        IS_ROTAMERIC,

        /**
         * Represents an amino acid that contains a terminal sp2-sp3 torsion.
         */
        NON_ROTAMERIC,

        /**
         * Represents an amino acid that does not have sidechain torsions.
         */
        HAS_NO_ROTAMERS,

        /**
         * Represents an unusual amino acid for which there will be no sidechain torsion data
         * in the Dunbrack database.
         */
        SPECIAL;
    }

    /**
     * Checks if this is a kind of proline.
     * @return true if this is a kind of proline
     */
    public boolean isProline()
    {
        if ( this == LPRO || this == DPRO )
            return true;
        return false;
    }
}
