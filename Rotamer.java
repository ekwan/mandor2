import java.util.*;
import java.io.*;
import org.apache.commons.math3.geometry.euclidean.threed.*;
import org.jgrapht.*;
import org.jgrapht.graph.*;
import org.jgrapht.alg.*;
import com.google.common.collect.*;
import java.util.concurrent.*;

/**
 * Represents a sidechain at a particular position in a given peptide.
 * Note: if we start to do amino acids with two substituents instead of one, we'll run into a problem with
 * the way rotamer is defined, since a mutation from one amino acid to another won't preserve the backbone.
 */
public class Rotamer implements Immutable
{
    /** the atoms in the sidechain */
    public final List<Atom> atoms;

    /** the position in the sequence of the original peptide */
    public final int sequenceIndex;

    /** the chi angles this rotamer has */
    public final List<Double> chis;

    /** the description field copied from Residue */
    public final String description;

    /**
     * Constructor.  Just copies the fields in the constructor.
     */
    public Rotamer(List<Atom> atoms, int sequenceIndex, List<Double> chis, String description)
    {
        if ( atoms == null )
            throw new NullPointerException("null atoms not allowed");
        if ( atoms.size() == 0 )
            throw new IllegalArgumentException("every rotamer must have some atoms");
        this.atoms = ImmutableList.copyOf(atoms);
        if ( sequenceIndex < 0 )
            throw new IllegalArgumentException("no negative sequence indices are allowed");
        this.sequenceIndex = sequenceIndex;
        if ( chis == null )
            throw new NullPointerException("null chis not allowed -- use empty list instead");
        this.chis = ImmutableList.copyOf(chis);
        if ( description == null || description.length() == 0 )
            throw new NullPointerException("must have a description");
        this.description = description;
    }

    /**
     * Mutates the given peptide to the specified rotamer. This will automatically
     * perform the mutation at the correct location in the original peptide. If the identity of the residue in
     * startingPeptide is wrong, we fix it first.
     * @param startingPeptide the peptide to make the mutation on
     * @param rotamer the rotamer we want in the resulting peptide
     * @return the mutated peptide
     */
    public static Peptide reconstitute(Peptide startingPeptide, Rotamer rotamer)
    {
        // mutate residue identity if necessary
        Peptide peptide = startingPeptide;
        if ( rotamer.sequenceIndex > startingPeptide.contents.size() - 1 )
            throw new IllegalArgumentException("sequence index out of bounds");
        Residue residue = startingPeptide.sequence.get(rotamer.sequenceIndex);
        if ( ! rotamer.description.equals(residue.description) )
            {
                ProtoAminoAcid protoAminoAcid = ProtoAminoAcidDatabase.getTemplate(rotamer.description);
                peptide = SidechainMutator.mutateSidechain(peptide, residue, protoAminoAcid);
                residue = peptide.sequence.get(rotamer.sequenceIndex);
           }
        
        // adjust chis (even if the chis are identical, we adjust anyways)
        if ( rotamer.chis.size() != residue.chis.size() )
            throw new IllegalArgumentException("chi size mismatch for " + rotamer.description + " at index " + rotamer.sequenceIndex);
        if (rotamer.chis.size() > 0)
            peptide = RotamerMutator.setChis(peptide, residue, rotamer.chis);
        return peptide;
    }

    /**
     * Mutates the given peptide to the specified rotamers.  No checks.
     * @param startingPeptide the starting peptide
     * @param rotamerList the rotamers to use
     * @return the mutated peptide
     */
    public static Peptide reconstitute(Peptide startingPeptide, List<Rotamer> rotamerList)
    {
        Peptide returnPeptide = startingPeptide;
        Set<Integer> checkSet = new HashSet<>();
        for (Rotamer r : rotamerList)
            {
                if ( checkSet.contains(r.sequenceIndex) )
                    throw new IllegalArgumentException("duplicate sequence index");
                else
                    checkSet.add(r.sequenceIndex);
            }
        for (Rotamer r : rotamerList)
            returnPeptide = reconstitute(returnPeptide, r);
        return returnPeptide;
    }

    @Override
    public boolean equals(Object obj)
    {
        if ( obj == null )
            return false;
        if ( obj == this )
            return true;
        if ( !(obj instanceof Rotamer) )
            return false;

        Rotamer r = (Rotamer)obj;
        if ( sequenceIndex == r.sequenceIndex &&
             chis.equals(r.chis) &&
             description.equals(r.description) )
            return true;
        return false;
    }

    @Override
    public int hashCode()
    {
        return Objects.hash(sequenceIndex, chis, description);
    }

    @Override
    public String toString()
    {
        String chiString = "";
        for (Double d : chis)
            chiString += String.format("%.1f, ", d);
        if ( chiString.length() > 0 )
            chiString = chiString.substring(0, chiString.length()-2);
        return String.format("[%d] %s: [%s]", sequenceIndex, description, chiString);
    }
}
