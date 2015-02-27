import java.util.*;
import java.io.*;
import org.apache.commons.math3.geometry.euclidean.threed.*;
import org.jgrapht.*;
import org.jgrapht.graph.*;
import org.jgrapht.alg.*;
import com.google.common.collect.*;
import java.util.concurrent.*;
import java.util.concurrent.atomic.*;

/**
 * This represents the possible rotamers at each position of a ground state peptide.
 */
public class VariableSequenceRotamerSpace extends RotamerSpace
{
    /**
     * Creates a VariableSequenceRotamerSpace just by copying fields.
     * @param peptide the peptide to pack the rotamers of
     * @param rotamerSpace the allowed rotamers at each position (an empty list means the position is fixed)
     * @param incompatiblePairs pairs of incompatible rotamers
     */
    public VariableSequenceRotamerSpace(Peptide peptide, List<List<Rotamer>> rotamerSpace, Set<RotamerPair> incompatiblePairs)
    {
        super(peptide, rotamerSpace, incompatiblePairs);
    }

    public VariableSequenceRotamerSpace(Peptide peptide, boolean includeHN)
    {
        super(peptide, includeHN, true);
    }

    public VariableSequenceRotamerSpace(Peptide peptide, boolean includeHN, boolean parallelize)
    {
        super(peptide, includeHN, parallelize);
    }

    /**
     * Figures out which rotamers are possible for a catalyst design.  Rotamers are placed at non-hairpin positions.
     * @param inputPeptide the peptide to analyze
     * @param includeHN whether to consider backbone HNs part of sidechains
     * @return the rotamer space (empty inner lists mean no variation is desired)
     */
    @Override
    public List<List<Rotamer>> getRotamerSpace(Peptide inputPeptide, boolean includeHN)
    {
        List<List<Rotamer>> rotamerSpace = new ArrayList<>(inputPeptide.sequence.size());
        for (Residue residue : inputPeptide.sequence)
            {
                if ( residue.isHairpin || residue.aminoAcid == AminoAcid.TS ||
                     residue.aminoAcid == AminoAcid.HIS || residue.aminoAcid == AminoAcid.ARG )
                    {
                        // place blank lists at positions with no rotamers
                        rotamerSpace.add(new ArrayList<Rotamer>());
                        //System.out.printf("Placed 0 rotamers at position %d (%s)\n", i, residue.description);
                    }
                else
                    {
                        // include only rotamers for the current residue
                        List<Rotamer> thisRotamers = RotamerFactory.generateRotamers(inputPeptide, residue, includeHN, null);
                        if ( thisRotamers.size() == 0 )
                            throw new IllegalArgumentException("expected at least one rotamer for this position!");
                        rotamerSpace.add(thisRotamers);
                        //System.out.printf("Placed %d rotamers at position %d (%s)\n", thisRotamers.size(), i, residue.description);
                    }
            }
        return rotamerSpace;
    }
}
