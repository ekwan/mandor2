import java.util.*;
import java.io.*;
import com.google.common.collect.*;
import java.util.concurrent.*;

/**
 * This singleton holds all of the ProtoAminoAcids.
 * Implemented as two parallel lists.
 */
public final class ProtoAminoAcidDatabase implements Singleton
{
    /**
     * Contains a list of all the possible amino acids.
     */
    public static final ImmutableList<AminoAcid> KEYS;

    /**
     * Contains all of the amino acid templates.
     * List indices parallel that of KEYS. 
     */
    public static final ImmutableList<ImmutableList<ProtoAminoAcid>> VALUES;

    // static initializer
    static
        {
            // make temporary lists
            List<AminoAcid> tempKeys = new LinkedList<>();
            List<List<ProtoAminoAcid>> tempValues = new LinkedList<>();
            Set<Atom> allAtoms = new HashSet<Atom>();

            // populate keys
            for (AminoAcid a : AminoAcid.values())
                {
                    tempKeys.add(a);
                    tempValues.add(new LinkedList<ProtoAminoAcid>());
                }

            // read all data from the library directory
            String directory = Settings.PROTOAMINOACID_DIRECTORY;
            int counter = 0;
            for (File f : new File(directory).listFiles())
                {
                    // parse all txt files
                    String filename = f.getName();
                    if ( filename.endsWith(".txt") && ! f.isDirectory() )
                        {
                            try
                                {
                                    // read file
                                    ProtoAminoAcidFile m = new ProtoAminoAcidFile(directory + filename);
                                    ProtoAminoAcid p = new ProtoAminoAcid(m);
                                    counter++;
                                    p = p.shift(counter);

                                    // check for duplicate atoms
                                    for (Atom a : p.molecule.contents)
                                        {
                                            if ( ! allAtoms.contains(a) )
                                                allAtoms.add(a);
                                            else
                                                throw new IllegalArgumentException("Duplicate atom: " + a.toString());
                                        }
                                    
                                    // add to database
                                    AminoAcid aminoAcid = p.residue.aminoAcid;
                                    int index = tempKeys.indexOf(aminoAcid);
                                    List<ProtoAminoAcid> thisList = tempValues.get(index);
                                    thisList.add(p);
                                }
                            catch (Exception e)
                                {
                                    System.out.println("Warning: error reading ProtoAminoAcid from: " + filename);
                                    e.printStackTrace();
                                }
                        }
                }

            // create immutable lists
            KEYS = ImmutableList.copyOf(tempKeys);
            
            List<ImmutableList<ProtoAminoAcid>> tempValues2 = new LinkedList<>();
            for (List<ProtoAminoAcid> l : tempValues)
                tempValues2.add( ImmutableList.copyOf(l) );
            VALUES = ImmutableList.copyOf(tempValues2);
        }

    /**
     * This class is not instantiable.
     */
    private ProtoAminoAcidDatabase()
    {
        throw new IllegalArgumentException("not instantiable");
    }

    /**
     * Given an amino acid, returns a list of all the ProtoAminoAcid templates
     * available for it.  If no templates are available, an empty list will be
     * returned.
     * @param aminoAcid the amino acid whose templates are desired
     * @return the immutable list of possible templates
     */
    public static ImmutableList<ProtoAminoAcid> getProtoAminoAcids(AminoAcid aminoAcid)
    {
        int index = KEYS.indexOf(aminoAcid);
        return VALUES.get(index);
    }

    /**
     * Returns a listing of which amino acids are available.
     * @return the description
     */
    public static String getDescription()
    {
        int total = 0;
        List<AminoAcid> empty = new LinkedList<>();
        String returnString = "== ProtoAminoAcidDatabase ==\n";
        
        // enumerate all possibilities
        for (int i=0; i < KEYS.size(); i++)
            {
                AminoAcid aminoAcid = KEYS.get(i);
                List<ProtoAminoAcid> currentList = VALUES.get(i);
                total += currentList.size();
                if ( currentList.size() == 1 )
                    {
                        returnString = returnString + currentList.get(0).residue.description;
                    }
                else if ( currentList.size() > 1 )
                    {
                        returnString = returnString + aminoAcid.fullName + ":\n";
                        for (int j=0; j < currentList.size(); j++)
                            {
                                ProtoAminoAcid p = currentList.get(j);
                                returnString = returnString + String.format("   %s", p.residue.description);
                                if ( j < currentList.size() - 1 )
                                    returnString = returnString + "\n";
                            }
                    }
                else
                    empty.add(aminoAcid);
                if ( currentList.size() > 0 )
                    returnString = returnString + "\n";
            }

        // list empty amino acids
        returnString = returnString + "\nAminoAcids without templates:\n";
        for (int i=0 ; i < empty.size(); i++)
            {
                AminoAcid a = empty.get(i);
                returnString = returnString + a.fullName;
                if ( i < empty.size() - 1 )
                    returnString = returnString + ", ";
            }
        returnString = returnString + "\n\n>>> " + total + " total entries";
        return returnString;
    }

    /**
     * Returns the ProtoAminoAcids whose descriptions match the given arguments.
     * If multiple ProtoAminoAcids match a description or no match is found,
     * an error will occur.  The comparison is case-insensitive.  The whole
     * description need not be used.  For example, "arg" would be enough to find
     * "standard_arginine," assuming there are no other descriptions containing "arg."
     * @param sequence a list of Strings containing the description of the requested amino acids
     * @return the list of ProtoAminoAcids in the same order as the input sequence
     */
    public static List<ProtoAminoAcid> getSpecificSequence(List<String> sequence)
    {
        // create a map of descriptions to ProtoAminoAcids
        Map<String,ProtoAminoAcid> masterMap = new HashMap<String,ProtoAminoAcid>();
        for (List<ProtoAminoAcid> list : VALUES)
            {
                for (ProtoAminoAcid p : list)
                    {
                        String description = p.residue.description.toLowerCase();
                        if ( masterMap.containsKey(description) )
                            throw new IllegalArgumentException("duplicate description found: " + description);
                        masterMap.put(description,p);
                    }
            }

        List<ProtoAminoAcid> returnList = new LinkedList<>();
        for (String s : sequence)
            {
                String thisDescription = s.toLowerCase();
                List<String> matches = new LinkedList<>();
                for (String key : masterMap.keySet())
                    {
                        if (key.indexOf(thisDescription) > -1)
                            matches.add(key);
                    }
                if ( matches.size() != 1 )
                    {
                        System.out.println(matches);
                        throw new IllegalArgumentException("invalid number of matches (" + matches.size() + " for " + s + ") -- cannot continue");
                    }
                ProtoAminoAcid hit = masterMap.get(matches.get(0));
                returnList.add(hit.shift(returnList.size()));
            }
        return ImmutableList.copyOf(returnList);
    }

    /** convenience method, hairpin set to true */
    public static List<ProtoAminoAcid> getSequence(String... args)
    {
        return getSequence(args, true);
    }

    /**
     * Creates a sequence from string input.  If there are multiple ProtoAminoAcids that
     * correspond to one amino acids, then one will be randomly selected.
     * If hairpin is turned off, glycine will explicitly not be set to a hairpin version.
     * @param args sequence of Strings representing the amino acid names (e.g., "pro", "gly")
     * @param hairpin true if you want auto-selection of the hairpin DPRO and GLY
     * @return the desired list of ProtoAminoAcids
     */
    public static List<ProtoAminoAcid> getSequence(String[] args, boolean hairpin)
    {
        List<ProtoAminoAcid> returnList = new LinkedList<>();
        for (int i=0; i < args.length; i++)
            {
                String s = args[i];
                AminoAcid aminoAcid = AminoAcid.getAminoAcid(s);
                List<ProtoAminoAcid> protoAminoAcidList = getProtoAminoAcids(aminoAcid);
                ProtoAminoAcid protoAminoAcid = null;
                if ( protoAminoAcidList.size() == 0 )
                    throw new IllegalArgumentException("no ProtoAminoAcids found for " + s);
                else if ( protoAminoAcidList.size() == 1 )
                    protoAminoAcid = protoAminoAcidList.get(0);
                else if ( protoAminoAcidList.size() > 1 )
                    {
                        // auto-select hairpin version following DPRO
                        if ( hairpin && i>0 && AminoAcid.getAminoAcid(args[i-1]) == AminoAcid.DPRO )
                            {
                                for (ProtoAminoAcid p : protoAminoAcidList)
                                    {
                                        if ( p.residue.isHairpin )
                                            {
                                                protoAminoAcid = p;
                                                //System.out.println("Auto-selected the hairpin version for " + s + ".");
                                            }
                                    }
                                if ( protoAminoAcid == null )
                                    {
                                        protoAminoAcid = protoAminoAcidList.get(0);
                                        System.out.println("Warning: No hairpin version found for " + s + "!!!");
                                    }
                            }
                        else
                            {
                                // tries to avoid hairpin versions, but 
                                // will select a hairpin version if nothing else is available
                                List<ProtoAminoAcid> pool = new LinkedList<>();
                                for (ProtoAminoAcid p : protoAminoAcidList)
                                    if ( !p.residue.isHairpin && p.residue.aminoAcid != AminoAcid.TS )
                                        pool.add(p);

                                ThreadLocalRandom random = ThreadLocalRandom.current();
                                if ( pool.size() == 0 )
                                    protoAminoAcid = protoAminoAcidList.get( random.nextInt(protoAminoAcidList.size()) );
                                else
                                    protoAminoAcid = pool.get( random.nextInt(pool.size()) );
                            }
                    }
                //System.out.println("Selected " + protoAminoAcid.residue.description);
                returnList.add(protoAminoAcid.shift(i));
            }
        return returnList;
    }

    /** for testing */
    public static void main(String[] args)
    {
        System.out.println(getDescription());
        for ( AminoAcid a : AminoAcid.values() )
        {   
            System.out.println();
            System.out.println(a);
            System.out.println(getProtoAminoAcids(a));
            System.out.println("\n\n");
        }
    }
}
