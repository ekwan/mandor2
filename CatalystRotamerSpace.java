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
 * This represents the possible rotamers at each position of a catalyst design.
 */
public class CatalystRotamerSpace extends RotamerSpace
{
    /** The possible mutation outcomes for non-interesting positions. */
    public static final List<ProtoAminoAcid> MUTATION_OUTCOMES;

    /** Static initiaizer. */
    static
    {
        // populate mutation outcomes
        List<ProtoAminoAcid> tempList = new ArrayList<>();
        List<AminoAcid> aminoAcids = ImmutableList.of(AminoAcid.ALA, AminoAcid.GLY, AminoAcid.VAL, AminoAcid.LEU, AminoAcid.ILE,
                                                      AminoAcid.PHE, AminoAcid.TYR, AminoAcid.TRP, AminoAcid.SER, AminoAcid.THR,
                                                      AminoAcid.ASN, AminoAcid.GLN, AminoAcid.ASP, AminoAcid.GLU);
        for (AminoAcid a : aminoAcids)
            tempList.addAll(ProtoAminoAcidDatabase.getProtoAminoAcids(a));
        MUTATION_OUTCOMES = ImmutableList.copyOf(tempList);
    }

    /**
     * Creates a CatalystRotamerSpace just by copying fields.
     * @param peptide the peptide to pack the rotamers of
     * @param rotamerSpace the allowed rotamers at each position (an empty list means the position is fixed)
     * @param incompatiblePairs pairs of incompatible rotamers
     */
    public CatalystRotamerSpace(Peptide peptide, List<List<Rotamer>> rotamerSpace, Set<RotamerPair> incompatiblePairs)
    {
        super(peptide, rotamerSpace, incompatiblePairs);
    }

    public CatalystRotamerSpace(Peptide peptide, boolean includeHN)
    {
        super(peptide, includeHN);
    }

    /**
     * Figures out which rotamers are possible for a catalyst design.  Rotamers are placed at non-hairpin
     * positions that do not have glycine.
     * @param inputPeptide the peptide to analyze
     * @param includeHN whether to consider backbone HNs part of sidechains
     * @return the rotamer space (empty inner lists mean no variation is desired)
     */
    @Override
    public List<List<Rotamer>> getRotamerSpace(Peptide inputPeptide, boolean includeHN)
    {
        int sequenceLength = inputPeptide.sequence.size();
        int forbiddenIndex = (sequenceLength/2) - 1;
        List<List<Rotamer>> rotamerSpace = new ArrayList<>(sequenceLength);
        for (int i=0; i < sequenceLength; i++)
            {
                Residue residue = inputPeptide.sequence.get(i);

                if ( i == forbiddenIndex || i == forbiddenIndex + 1 || residue.aminoAcid != AminoAcid.GLY )
                    {
                        // place blank lists at positions with no rotamers
                        rotamerSpace.add(new ArrayList<Rotamer>());
                    }
                else
                    {
                        // mutate to all possible amino acids and place rotamers for each
                        List<Rotamer> allRotamers = new ArrayList<>();
                        for (ProtoAminoAcid protoAminoAcid : MUTATION_OUTCOMES)
                            {
                                Peptide mutatedPeptide = SidechainMutator.mutateSidechain(inputPeptide, residue, protoAminoAcid);
                                Residue mutatedResidue = mutatedPeptide.sequence.get(i);
                                List<Rotamer> thisRotamers = RotamerFactory.generateRotamers(mutatedPeptide, mutatedResidue, includeHN, null);
                                allRotamers.addAll(thisRotamers);
                            }
                        if ( allRotamers.size() == 0 )
                            throw new IllegalArgumentException("expected at least one rotamer for this position!");
                        rotamerSpace.add(allRotamers);
                    }
            }
        return rotamerSpace;
    }

    public static void main(String[] args)
    {
        DatabaseLoader.go();
        List<Peptide> sheets = BetaSheetGenerator.generateSheets(6, 5, 10000, 0.01);
        Peptide peptide = sheets.get(0);
        CatalystRotamerSpace catalystRotamerSpace = new CatalystRotamerSpace(peptide, true);
    }
}
