import java.util.*;
import com.google.common.collect.*;
import org.jgrapht.*;
import org.jgrapht.graph.*;

/**
 * This class converts the atom types in a Peptide from one that
 * has the histidine and hydroxyl far apart to one that has them
 * close together.  This involves adding a bond to the connectivity table.
 * We can also undo this change.
 */
public class HydrogenBondMutator implements Singleton, Mutator
{
    public static final List<Pair<Integer,Integer>> HID_SEPARATED_ATOM_TYPES;
    public static final List<Pair<Integer,Integer>> HID_CLOSE_ATOM_TYPES;
    public static final List<Pair<Integer,Integer>> HIE_SEPARATED_ATOM_TYPES;
    public static final List<Pair<Integer,Integer>> HIE_CLOSE_ATOM_TYPES;
    
    /**
     * the new hydrogen bond will be formed from this close contact atom type (AMOEBA)
     */
    public static final int HBOND_FROM_ATOM_TYPE = 416;

    /**
     * the new hydrogen bond will be formed to this close contact atom type (AMOEBA)
     */
    public static final int HBOND_TO_ATOM_TYPE = 421;

    static
    {
        // check that the hydrogen bond is being formed between atoms of different types
        if ( HBOND_FROM_ATOM_TYPE == HBOND_TO_ATOM_TYPE )
            throw new IllegalArgumentException("must form hydrogen bond between atoms of different types");

        // initialize parallel lists
        List<Pair<Integer,Integer>> tempCommonSeparatedList = new LinkedList<>();
        List<Pair<Integer,Integer>> tempCommonCloseList = new LinkedList<>();
        List<Pair<Integer,Integer>> tempHIDSeparatedList = new LinkedList<>();
        List<Pair<Integer,Integer>> tempHIDCloseList = new LinkedList<>();
        List<Pair<Integer,Integer>> tempHIESeparatedList = new LinkedList<>();
        List<Pair<Integer,Integer>> tempHIECloseList = new LinkedList<>();

        // regular transition state atom type conversion
        // these entries are common to both HID and HIE
        // first entry in pair is the AMOEBA type, second is the OPLS type
        tempCommonSeparatedList.add(new Pair<Integer,Integer>(401,401)); tempCommonCloseList.add(new Pair<Integer,Integer>(421,421));
        tempCommonSeparatedList.add(new Pair<Integer,Integer>(402,402)); tempCommonCloseList.add(new Pair<Integer,Integer>(422,422));
        tempCommonSeparatedList.add(new Pair<Integer,Integer>(403,403)); tempCommonCloseList.add(new Pair<Integer,Integer>(423,423));
        tempCommonSeparatedList.add(new Pair<Integer,Integer>(404,404)); tempCommonCloseList.add(new Pair<Integer,Integer>(424,424));
        tempCommonSeparatedList.add(new Pair<Integer,Integer>(405,405)); tempCommonCloseList.add(new Pair<Integer,Integer>(425,425));
        tempCommonSeparatedList.add(new Pair<Integer,Integer>(406,406)); tempCommonCloseList.add(new Pair<Integer,Integer>(426,426));
        tempCommonSeparatedList.add(new Pair<Integer,Integer>(407,407)); tempCommonCloseList.add(new Pair<Integer,Integer>(427,427));
        tempCommonSeparatedList.add(new Pair<Integer,Integer>(408,408)); tempCommonCloseList.add(new Pair<Integer,Integer>(428,428));
        tempCommonSeparatedList.add(new Pair<Integer,Integer>(409,409)); tempCommonCloseList.add(new Pair<Integer,Integer>(429,429));
        tempCommonSeparatedList.add(new Pair<Integer,Integer>(410,410)); tempCommonCloseList.add(new Pair<Integer,Integer>(430,430));
        tempCommonSeparatedList.add(new Pair<Integer,Integer>(411,411)); tempCommonCloseList.add(new Pair<Integer,Integer>(431,431));

        tempHIDSeparatedList.addAll(tempCommonSeparatedList);
        tempHIDCloseList.addAll(tempCommonCloseList);
        tempHIESeparatedList.addAll(tempCommonSeparatedList);
        tempHIECloseList.addAll(tempCommonCloseList);

        // HID histidine to close contact atom types
        tempHIDSeparatedList.add(new Pair<Integer,Integer>(119,153)); tempHIDCloseList.add(new Pair<Integer,Integer>(414,414));
        tempHIDSeparatedList.add(new Pair<Integer,Integer>(122,152)); tempHIDCloseList.add(new Pair<Integer,Integer>(412,412));
        tempHIDSeparatedList.add(new Pair<Integer,Integer>(123, 12)); tempHIDCloseList.add(new Pair<Integer,Integer>(413,413));
        tempHIDSeparatedList.add(new Pair<Integer,Integer>(126,156)); tempHIDCloseList.add(new Pair<Integer,Integer>(416,416));
        tempHIDSeparatedList.add(new Pair<Integer,Integer>(124,151)); tempHIDCloseList.add(new Pair<Integer,Integer>(419,419));
        tempHIDSeparatedList.add(new Pair<Integer,Integer>(125, 12)); tempHIDCloseList.add(new Pair<Integer,Integer>(420,420));
        tempHIDSeparatedList.add(new Pair<Integer,Integer>(120,148)); tempHIDCloseList.add(new Pair<Integer,Integer>(417,417));
        tempHIDSeparatedList.add(new Pair<Integer,Integer>(121,149)); tempHIDCloseList.add(new Pair<Integer,Integer>(418,418));

        // HIE histidine to close contact atom types
        tempHIESeparatedList.add(new Pair<Integer,Integer>(129,152)); tempHIECloseList.add(new Pair<Integer,Integer>(412,412));
        tempHIESeparatedList.add(new Pair<Integer,Integer>(131,153)); tempHIECloseList.add(new Pair<Integer,Integer>(414,414));
        tempHIESeparatedList.add(new Pair<Integer,Integer>(132, 12)); tempHIECloseList.add(new Pair<Integer,Integer>(415,415));
        tempHIESeparatedList.add(new Pair<Integer,Integer>(135,148)); tempHIECloseList.add(new Pair<Integer,Integer>(417,417));
        tempHIESeparatedList.add(new Pair<Integer,Integer>(136,149)); tempHIECloseList.add(new Pair<Integer,Integer>(418,418));
        tempHIESeparatedList.add(new Pair<Integer,Integer>(133,151)); tempHIECloseList.add(new Pair<Integer,Integer>(419,419));
        tempHIESeparatedList.add(new Pair<Integer,Integer>(134, 12)); tempHIECloseList.add(new Pair<Integer,Integer>(420,420));
        tempHIESeparatedList.add(new Pair<Integer,Integer>(130,156)); tempHIECloseList.add(new Pair<Integer,Integer>(416,416));

        // make immutable copies
        HID_SEPARATED_ATOM_TYPES  = ImmutableList.copyOf(tempHIDSeparatedList);
        HID_CLOSE_ATOM_TYPES      = ImmutableList.copyOf(tempHIDCloseList);
        HIE_SEPARATED_ATOM_TYPES  = ImmutableList.copyOf(tempHIESeparatedList);
        HIE_CLOSE_ATOM_TYPES      = ImmutableList.copyOf(tempHIECloseList);
        
        // check list sizes
        if ( HID_SEPARATED_ATOM_TYPES.size() != HID_CLOSE_ATOM_TYPES.size() )
            throw new IllegalArgumentException(String.format("HID list size mismatch %d %d", HID_SEPARATED_ATOM_TYPES.size(), HID_CLOSE_ATOM_TYPES.size()));
        if ( HIE_SEPARATED_ATOM_TYPES.size() != HIE_CLOSE_ATOM_TYPES.size() )
            throw new IllegalArgumentException(String.format("HID list size mismatch %d %d", HIE_SEPARATED_ATOM_TYPES.size(), HIE_CLOSE_ATOM_TYPES.size()));

        // check every key has a value and there are no duplicate keys
        // (duplicate values might occur)
        for (int i=0; i < HID_SEPARATED_ATOM_TYPES.size(); i++)
            {
                Pair<Integer,Integer> from = HID_SEPARATED_ATOM_TYPES.get(i);
                Pair<Integer,Integer> to   = HID_CLOSE_ATOM_TYPES.get(i);
                if ( from == null || to == null ||
                     from.getFirst() < 1 || from.getSecond() < 1 ||
                     to.getFirst() < 1   || to.getSecond() < 1 ||
                     from.equals(to) )
                    throw new IllegalArgumentException("invalid entry in HID table");
                int occurrences = Collections.frequency(HID_SEPARATED_ATOM_TYPES, from);
                if ( occurrences > 1 )
                    throw new IllegalArgumentException(String.format("duplicate HID key found (%d instances, from %s to %s", occurrences, from.toString(), to.toString()));
            }

        for (int i=0; i < HIE_SEPARATED_ATOM_TYPES.size(); i++)
            {
                Pair<Integer,Integer> from = HIE_SEPARATED_ATOM_TYPES.get(i);
                Pair<Integer,Integer> to   = HIE_CLOSE_ATOM_TYPES.get(i);
                if ( from == null || to == null ||
                     from.getFirst() < 1 || from.getSecond() < 1 ||
                     to.getFirst() < 1   || to.getSecond() < 1 ||
                     from.equals(to) )
                    throw new IllegalArgumentException("invalid entry in HIE table");
                int occurrences = Collections.frequency(HIE_SEPARATED_ATOM_TYPES, from);
                if ( occurrences > 1 )
                    throw new IllegalArgumentException(String.format("duplicate HIE key found (%d instances, from %s to %s", occurrences, from.toString(), to.toString()));
            }
    }

    /** class is not instantiable */
    private HydrogenBondMutator()
    {
        throw new IllegalArgumentException("not instantiable");
    }

    /**
     * Returns the separated and close atom types appropriate for this peptide.
     * @param peptide the peptide
     * @return first is separated second is close
     */
    public static Pair< List<Pair<Integer,Integer>>, List<Pair<Integer,Integer>> > getAtomTypeLists(Peptide peptide)
    {
        // determine whether we are dealing with HID or HIE
        String description = null;
        for ( Residue r : peptide.sequence )
            {
                if ( r.description.indexOf("histidine") > -1 )
                    {
                        if ( description == null )
                            description = r.description;
                        else
                            throw new IllegalArgumentException("cannot handle two histidines");
                    }
            }
        if ( description == null )
            throw new NullPointerException("must be exactly one histidine in the molecule to do a mutation");

        if ( description.indexOf("histidine_hd") > -1 )
            return new Pair<>(HID_SEPARATED_ATOM_TYPES, HID_CLOSE_ATOM_TYPES);
        else if ( description.indexOf("histidine_he") > -1 )
            return new Pair<>(HIE_SEPARATED_ATOM_TYPES, HIE_CLOSE_ATOM_TYPES);
        else
            throw new IllegalArgumentException("unknown histidine type");
    }

    /**
     * Converts the atom types in a peptide from ones where the histidine and
     * hydroxyl are far apart to ones where they are close together.  Also
     * adds the new hydrogen bond.
     * @param peptide the original peptide
     * @return the new peptide
     */
    public static Peptide mutate(Peptide peptide)
    {
        // adjust atom types creating new atom map
        Pair< List<Pair<Integer,Integer>>, List<Pair<Integer,Integer>> > listPair = getAtomTypeLists(peptide);
        List<Pair<Integer,Integer>> separatedTypes = listPair.getFirst();
        List<Pair<Integer,Integer>> closeTypes = listPair.getSecond();
        HashMap<Atom,Atom> atomMap = new HashMap<>();
        for (Atom a : peptide.contents)
            {
                Pair<Integer,Integer> currentType = new Pair<>(a.type1, a.type2);
                if ( separatedTypes.contains(currentType) )
                    {
                        int index = separatedTypes.indexOf(currentType);
                        Pair<Integer,Integer> targetType = closeTypes.get(index);
                        Atom newAtom = a.setAtomType(targetType.getFirst(), targetType.getSecond());
                        atomMap.put(a, newAtom);
                    }
            }
        Peptide adjustedPeptide = peptide.moveAtoms2(atomMap);

        // locate the atoms to form a hydrogen bond between
        Atom fromAtom = null;
        Atom toAtom = null;
        for (Atom a : adjustedPeptide.contents)
            {
                if ( a.type1 == HBOND_FROM_ATOM_TYPE )
                    {
                        if ( fromAtom == null )
                            fromAtom = a;
                        else
                            throw new IllegalArgumentException("error, duplicate from atom type when forming hydrogen bond!");
                    }

                if ( a.type1 == HBOND_TO_ATOM_TYPE )
                    {
                        if ( toAtom == null )
                            toAtom = a;
                        else
                            throw new IllegalArgumentException("error, duplicate to atom type when forming hydrogen bond!");
                    }
            }
        if ( fromAtom == null || toAtom == null )
            throw new NullPointerException("couldn't find the atoms to make a hydrogen bond between");

        // add the new bond
        if ( adjustedPeptide.connectivity.containsEdge(fromAtom,toAtom) )
            throw new IllegalArgumentException("target hydrogen bond already exists");
        DefaultWeightedEdge newEdge = adjustedPeptide.connectivity.addEdge(fromAtom,toAtom);
        adjustedPeptide.connectivity.setEdgeWeight(newEdge, 1.0);

        return adjustedPeptide;
    }

    /**
     * Converts the atom types in a peptide from the close contact types to the
     * separated types.  Also removes the new hydrogen bond.
     * @param peptide the original peptide
     * @return the new peptide
     */
    public static Peptide unmutate(Peptide peptide)
    {
        // locate the atoms to remove a hydrogen bond between
        Atom fromAtom = null;
        Atom toAtom = null;
        for (Atom a : peptide.contents)
            {
                if ( a.type1 == HBOND_FROM_ATOM_TYPE )
                    {
                        if ( fromAtom == null )
                            fromAtom = a;
                        else
                            throw new IllegalArgumentException("error, duplicate from atom type when deleting hydrogen bond!");
                    }

                if ( a.type2 == HBOND_TO_ATOM_TYPE )
                    {
                        if ( toAtom == null )
                            toAtom = a;
                        else
                            throw new IllegalArgumentException("error, duplicate to atom type when deleting hydrogen bond!");
                    }
            }
        if ( fromAtom == null || toAtom == null )
            throw new NullPointerException("couldn't find the atoms to delete a hydrogen bond with");

        // delete the hydrogen bond
        Peptide adjustedPeptide = peptide.moveAtoms2(new HashMap<Atom,Atom>());
        DefaultWeightedEdge edge = adjustedPeptide.connectivity.getEdge(fromAtom,toAtom);
        if ( edge == null )
            throw new NullPointerException("found the hydrogen-bonded atoms, but they're not bonded");
        adjustedPeptide.connectivity.removeEdge(edge);

        // adjust atom types creating new atom map
        Pair< List<Pair<Integer,Integer>>, List<Pair<Integer,Integer>> > listPair = getAtomTypeLists(peptide);
        List<Pair<Integer,Integer>> separatedTypes = listPair.getFirst();
        List<Pair<Integer,Integer>> closeTypes = listPair.getSecond();
        
        HashMap<Atom,Atom> atomMap = new HashMap<>();
        for (Atom a : adjustedPeptide.contents)
            {
                Pair<Integer,Integer> currentType = new Pair<>(a.type1, a.type2);
                if ( closeTypes.contains(currentType) )
                    {
                        int index = closeTypes.indexOf(currentType);
                        Pair<Integer,Integer> targetType = separatedTypes.get(index);
                        Atom newAtom = a.setAtomType(targetType.getFirst(), targetType.getSecond());
                        atomMap.put(a, newAtom);
                    }
            }        
        adjustedPeptide = adjustedPeptide.moveAtoms2(atomMap);

        // transfer the EnergyBreakdown
        //EnergyBreakdown energyBreakdown = peptide.energyBreakdown;
        //adjustedPeptide = adjustedPeptide.setEnergyBreakdown(energyBreakdown);
        return adjustedPeptide;
    }

    /**
     * Mutate a whole bunch of peptides at once.
     * @param peptides the peptides to mutate
     * @param the mutated peptides
     */
    public static List<Peptide> mutate(List<Peptide> peptides)
    {
        List<Peptide> returnList = new ArrayList<>(peptides.size());
        for (Peptide p : peptides)
            returnList.add(mutate(p));
        return returnList;
    }

    /**
     * Mutate a whole bunch of peptides at once.
     * @param peptides the peptides to mutate
     * @param the mutated peptides
     */
    public static List<Peptide> unmutate(List<Peptide> peptides)
    {
        List<Peptide> returnList = new ArrayList<>(peptides.size());
        for (Peptide p : peptides)
            returnList.add(unmutate(p));
        return returnList;
    }
}
