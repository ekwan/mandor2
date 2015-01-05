
// This class represents an amino acid.

import java.util.*;
import java.io.*;
import com.google.common.collect.*;
import java.util.concurrent.*;

/**
 * Represents the standard amino acids with different categories
 * for cis and trans proline.  Backbone-dependent rotamer data is read
 * from the Dunbrack library.  We use traditional rotamers for amino acids
 * containing all sp3-sp3 bonds in their sidechains.  That is, we use ordered
 * tuples (X<sub>1</sub>, X<sub>2</sub>, ..., X<sub>n</sub>) to represent the
 * backbone torsion angles, where n is the number of sidechain torsions.  (OHs
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
    /** all amino acids are L-configured unless otherwise noted */
    ALA("Ala",             "alanine",                      RotamerType.HAS_NO_ROTAMERS, Chirality.L),
    GLY("Gly",             "glycine",                      RotamerType.HAS_NO_ROTAMERS, Chirality.ACHIRAL),
    VAL("Val",             "valine",                       RotamerType.IS_ROTAMERIC, Chirality.L),
    LEU("Leu",             "leucine",                      RotamerType.IS_ROTAMERIC, Chirality.L),
    ILE("Ile",             "isoleucine",                   RotamerType.IS_ROTAMERIC, Chirality.L),

    LPRO("L-Pro",          "L-proline",                    RotamerType.IS_ROTAMERIC, Chirality.L),
    DPRO("D-Pro",          "D-proline",                    RotamerType.SPECIAL, Chirality.D),

    PHE("Phe",             "phenylalanine",                RotamerType.NON_ROTAMERIC, Chirality.L), 
    TYR("Tyr",             "tyrosine",                     RotamerType.NON_ROTAMERIC, Chirality.L),
    TRP("Trp",             "tryptophan",                   RotamerType.NON_ROTAMERIC, Chirality.L),
    SER("Ser",             "serine",                       RotamerType.IS_ROTAMERIC, Chirality.L),
    THR("Thr",             "threonine",                    RotamerType.IS_ROTAMERIC, Chirality.L),
    CYS("Cys",             "cysteine",                     RotamerType.IS_ROTAMERIC, Chirality.L),
    MET("Met",             "methionine",                   RotamerType.IS_ROTAMERIC, Chirality.L),
    ASN("Asn",             "aspargine",                    RotamerType.NON_ROTAMERIC, Chirality.L),
    GLN("Gln",             "glutamine",                    RotamerType.NON_ROTAMERIC, Chirality.L),
    LYS("Lys",             "lysine",                       RotamerType.IS_ROTAMERIC, Chirality.L),
    ARG("Arg",             "arginine",                     RotamerType.IS_ROTAMERIC, Chirality.L),
    HIS("His",             "histidine",                    RotamerType.NON_ROTAMERIC, Chirality.L),
    ASP("Asp",             "aspartate",                    RotamerType.NON_ROTAMERIC, Chirality.L),
    GLU("Glu",             "glutamate",                    RotamerType.NON_ROTAMERIC, Chirality.L);
    
    // fields

    /**
     * An abbreviation like "Ala".
     */
    public final String shortName;

    /**
     * A full name like "alanine".
     */
    public final String fullName;

    /**
     * Indicates the kinds of sidechain torsions present.
     */
    public final RotamerType rotamerType;

    /**
     * The filename containing the data for this residue.
     */
    public final String filename;

    /** 
     * The chirality of an amino acid.
     */
    public final Chirality chirality;

    /**
     * Contains the rotamer library data for this amino acid.
     */
    private SidechainRotamerLibrary library;

    // enum constructor
    AminoAcid(String shortName, String fullName, RotamerType rotamerType, Chirality chirality)
    {
        this.fullName = fullName;
        this.shortName = shortName;
        this.rotamerType = rotamerType;
        this.chirality = chirality;

        // determine filename
        if ( rotamerType == RotamerType.IS_ROTAMERIC )
	    {
		if (shortName.equals("L-Cpr"))
		    {
			filename = Settings.ROTAMER_LIBRARY_DIRECTORY + "cpr" + ".bbdep.rotamers.lib";
		    }
		else if (shortName.equals("L-Tpr"))
		    {
			filename = Settings.ROTAMER_LIBRARY_DIRECTORY + "tpr" + ".bbdep.rotamers.lib";
		    }
		else
		    {
			filename = Settings.ROTAMER_LIBRARY_DIRECTORY + shortName.toLowerCase() + ".bbdep.rotamers.lib";
		    }
	    }
        else if ( rotamerType == RotamerType.NON_ROTAMERIC )
	    filename = Settings.ROTAMER_LIBRARY_DIRECTORY + shortName.toLowerCase() + ".bbdep.densities.lib";
        else if ( rotamerType == RotamerType.HAS_NO_ROTAMERS || rotamerType == RotamerType.SPECIAL )
            filename = "";
        else
            throw new IllegalArgumentException("Unrecognized RotamerType in AminoAcid constructor!");

	    //call constructor of side chain rotamer library
    }

    /**
     * Method to extract amino acid type from string
	 * @param input that contains a string with an amino acid name
	 * @return the amino acid that corresponds to the string input
    */
    public static AminoAcid getAminoAcid(String input)
    {
        return valueOf(input.toUpperCase());
    }

    /**
     * Returns an amino acid that isn't special for testing purposes.
     * Possible return list does not include lysine, cysteine, and methionine.
     * @return the random amino acid
     */
    public static String getRandom()
    {
        List<AminoAcid> pool = new LinkedList<>();
        for (AminoAcid a : values())
            pool.add(a);

        pool.remove(AminoAcid.LYS);
        pool.remove(AminoAcid.CYS);
        pool.remove(AminoAcid.MET);

        pool.remove(AminoAcid.DPRO);
        pool.remove(AminoAcid.CPR);
        pool.remove(AminoAcid.TPR);

        ThreadLocalRandom random = ThreadLocalRandom.current();
        return pool.get(random.nextInt(pool.size())).name();
    }

    /**
     * Returns a brief description of this amino acid.
     * @return the textual description of this amino acid
     */
    @Override
    public String toString()
    {
        return shortName;
    }

    /**
     * Returns the filename associated with an amino acid
     * @return the String containing the filename
     */
    public String getFilename()
    {
	    return filename;
    }

    /**
     * Indicates whether the amino acid can be represented by standard rotamers,
     * is non-rotameric, or contains no rotatable bonds at all.  For a full description,
     * see the <a href="http://dunbrack.fccc.edu/bbdep2010/">Dunbrack backbone-dependent
     *   rotamer library</a> page.
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
         * Represents ALL of the amino acids for the Ramachandran library.
         */
        SPECIAL;
    }

    /**
     * Indicates the chirality of the amino acid. Options include L, D, and achiral (glycine).
     */
    public enum Chirality
    {
        /**
         * Represents a regular L amino acid
         */
        L,

        /**
         * Represents a D amino acid, which is most commonly D-proline
         */
        D,

        /**
         * Amino acid that lacks chirarlity. Glyicne is the only achiral amino acid. 
         */
        ACHIRAL;

    }

    /**
     * Gives a random set of sidechain torsion angles for this amino acid.
     * @param psi the backbone angle
     * @param phi the backbone angle
     * @return the torsion angles X1, X2, ... as an ordered list in degrees
     */
    
    public List<Double> getRandomRotamer(Double phi, Double psi)
    {
	    if (rotamerType == RotamerType.HAS_NO_ROTAMERS)
	        throw new IllegalArgumentException("Non rotameric amino acid");
	    else if (rotamerType == RotamerType.IS_ROTAMERIC)
            {
	            RotamericLibrary rotLib = new RotamericLibrary(this);
	            DiscreteProbabilityDistribution<List<Double>> dpd = rotLib.get(phi,psi);
	            return dpd.getRandom();
	        }
	    else if (rotamerType == RotamerType.NON_ROTAMERIC)
            {
	            NonRotamericLibrary nRotLib = new NonRotamericLibrary(this);
	            DiscreteProbabilityDistribution<NonRotamericLibrary.NonRotamericAngles> dpd1 = nRotLib.get(phi,psi);
	            NonRotamericLibrary.NonRotamericAngles nrA = dpd1.getRandom();
	            Double lastChi = nrA.getDPD().getRandom();
                List<Double> returnList = new LinkedList<>(nrA.getRotamericAngles());
                returnList.add(lastChi);
                return returnList;
	        }
        // should be unreachable
	    throw new IllegalArgumentException("cannot get a random rotamer for this kind of amino acid");
    }
    
    /**
     * Return what kind of rotamers this amino acid has.
     * @return the rotamer type
     */
    public RotamerType getRotamerType()
    {
	    return rotamerType;
    }

    /**
     * Returns true if this is a kind of proline.
     * @return true if proline
     */
    public boolean isProline()
    {
        if ( fullName.indexOf("proline") > -1 )
            return true;
        return false;
    }

    /**
     * Returns an amino acid compatible with the Ramachandran library.
     * @return the appropriate key for the database
     */
    public AminoAcid getCompatible()
    {
        if ( this == DPRO )
            return TPR;
        else
            return this;
    }

    // for testing
    public static void main(String[] args)
    {
        //System.out.println(AminoAcid.getTotalRotamers());

	    //System.out.println(AminoAcid.getRotamer(AminoAcid.Gln, 120.0, 120.0));
    }
}
