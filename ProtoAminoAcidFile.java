import java.io.*;
import java.util.*;
import com.google.common.collect.*;

/**
 * Reads files that contain information about special atoms, torsions, etc. for ProtoAminoAcids.
 * toString, equals, and hashCode not implemented.
 */
public class ProtoAminoAcidFile extends OutputFileFormat implements Serializable, Immutable
{
    /** For serialization. */
    public static final long serialVersionUID = 1L;

    /** Defined as the ordered pair (terminal acetyl amide carbon, terminal acetyl amide nitrogen) */
    public final Pair<Atom,Atom> NStickyConnection;

    /** Defined as the ordered pair (amide carbonyl carbon, NH2 nitrogen) */
    public final Pair<Atom,Atom> CStickyConnection;

    /** The template molecule. */
    public final Molecule molecule;

    /** The residue metadata.  Description is trimmed for whitespace and taken only in lowercase. */
    public final Residue residue;

    /** public constructor */
    public ProtoAminoAcidFile(String filename)
    {
        super(filename);

        // check that this is a valid file
        String headerString = stringRepresentation.split("\n")[0];
        if ( headerString.indexOf("Amino Acid Metadata") == -1 )
            throw new IllegalArgumentException("this does not appear to be an amino acid metadata file");

        // setup fields
        Pair<Atom,Atom> tempNStickyConnection = null;
        Pair<Atom,Atom> tempCStickyConnection = null;
        Molecule        tempMolecule          = null;
        Residue         tempResidue           = null;

        // parse molecule file
        // first we need the name of the file
        String templateFilename = null;
        for (List<String> line : fileContents)
            {
                if ( line.size() > 2 && line.get(0).toLowerCase().equals("xyz_file") )
                    {
                        if ( templateFilename != null )
                            throw new IllegalArgumentException("duplicate xyz_file field");
                        templateFilename = line.get(1);
                    }
            }
        if ( templateFilename == null )
            throw new NullPointerException("xyz filename not found");
        templateFilename = Settings.PROTOAMINOACID_DIRECTORY + templateFilename;
        TinkerXYZOutputFile templateFile = new TinkerXYZOutputFile(templateFilename);
        tempMolecule = templateFile.molecule;

        // setup some more temporary fields
        AminoAcid            aminoAcid           = null;
        String               description         = null;
        ProtoTorsion         phi                 = null;
        ProtoTorsion         psi                 = null;
        ProtoTorsion         omega               = null;
        List<String>         chiTitles           = new ArrayList<>();
        List<ProtoTorsion>   chis                = new ArrayList<>();
        Atom                 HN                  = null;
        Atom                 N                   = null;
        Atom                 O                   = null;
        Atom                 C                   = null;
        Atom                 CA                  = null;
        Atom                 HA                  = null;
        Pair<Atom,Atom>      prochiralConnection = null;
        List<Atom>           atoms               = new ArrayList<>(tempMolecule.contents);
        Boolean              isHairpin           = false;
        List<Integer>        OPLSatomTypes       = new ArrayList<>(tempMolecule.contents.size());
        List<Double>         surfaceTensions     = new ArrayList<>(tempMolecule.contents.size());

        // parse metadata file
        for (List<String> line : fileContents)
            {
                // check for blank or comment lines
                if ( line.size() < 2 || line.get(0).startsWith("!") )
                    continue;
                
                String firstField = line.get(0).toLowerCase();

                if ( firstField.equals("amino_acid") )
                    aminoAcid = AminoAcid.getAminoAcid(line.get(1));
                else if ( firstField.equals("description") )
                    {
                        description = "";
                        for ( int i=1; i < line.size(); i++ )
                            {
                                description += line.get(i);
                                if ( i < line.size() )
                                    description += " ";
                            }
                        description = description.toLowerCase().trim();
                    }
                else if ( firstField.equals("phi") )
                    {
                        if ( phi != null )
                            throw new IllegalArgumentException("phi already set");
                        phi = parseProtoTorsion(tempMolecule, line);
                   }
                else if ( firstField.equals("psi") )
                    {
                        if ( psi != null )
                            throw new IllegalArgumentException("psi already set");
                        psi = parseProtoTorsion(tempMolecule, line);
                    }
                else if ( firstField.equals("omega") )
                    {
                        if ( omega != null )
                            throw new IllegalArgumentException("omega already set");
                        omega = parseProtoTorsion(tempMolecule, line);
                    }
               else if ( firstField.startsWith("chi") )
                    {
                        if ( chiTitles.contains(firstField) )
                            throw new IllegalArgumentException("duplicate chi field: " + firstField);
                        
                        chiTitles.add(firstField);
                        ProtoTorsion thisTorsion = null;
                        if ( line.size() >= 6 && ( line.get(5).toUpperCase().equals("A") ) )
                            thisTorsion = parseAtomQuadruple(tempMolecule, line);
                        else
                            thisTorsion = parseProtoTorsion(tempMolecule, line);

                        if ( thisTorsion != null )
                            {
                                chis.add(thisTorsion);
                                int expectedSize = Integer.parseInt(firstField.replaceAll("[^0-9]",""));
                                if ( chis.size() != expectedSize )
                                    throw new IllegalArgumentException("are the chis out of order?");
                            }
                        else
                            throw new IllegalArgumentException("invalid chi field: " + line.toString());
                    }
                else if ( firstField.equals("atomhn") )
                    {
                        if ( HN != null )
                            throw new IllegalArgumentException("HN already set");
                        HN = parseSingleAtom(tempMolecule, line, ImmutableList.of(4, 10));
                    }

                // in the section below, the lists are of the autodetected amoeba types 
                else if ( firstField.equals("atomn") )
                    {
                        if ( N != null )
                            throw new IllegalArgumentException("N already set");
                        N = parseSingleAtom(tempMolecule, line, ImmutableList.of(1, 7, 50));
                    }
                else if ( firstField.equals("atomo") )
                    {
                        if ( O != null )
                            throw new IllegalArgumentException("O already set");
                        O = parseSingleAtom(tempMolecule, line, ImmutableList.of(5, 11, 53));
                    }
                else if ( firstField.equals("atomc") )
                    {
                        if ( C != null )
                            throw new IllegalArgumentException("C already set");
                        C = parseSingleAtom(tempMolecule, line, ImmutableList.of(3, 9, 52));
                    }
                else if ( firstField.equals("atomca") )
                    {
                        if ( CA != null )
                            throw new IllegalArgumentException("CA already set");
                        CA = parseSingleAtom(tempMolecule, line, ImmutableList.of(2, 8, 51));
                    }
                else if ( firstField.equals("atomha") )
                    {
                        if ( HA != null )
                            throw new IllegalArgumentException("HA already set");
                        HA = parseSingleAtom(tempMolecule, line, ImmutableList.of(6, 12, 54 ));
                    }

                // read the atom pairs
                else if ( firstField.equals("nstickyconnection") )
                    tempNStickyConnection = parsePair(tempMolecule, line);

                else if ( firstField.equals("cstickyconnection") )
                    tempCStickyConnection = parsePair(tempMolecule, line);

                else if ( firstField.equals("prochiral_connection") )
                    {
                        if ( line.size() >= 4 && line.get(3).toUpperCase().equals("A") )
                            prochiralConnection = parseAtomPair(tempMolecule, line);
                        else
                            prochiralConnection = parsePair(tempMolecule, line);
                    }
 
                // read the OPLS atom types
                else if ( firstField.equals("oplstype") )
                    {
                        if ( line.size() < 3 )
                            throw new IllegalArgumentException("unexpected number of fields for OPLS type\n" + line.toString());
                        int atomNumber = Integer.parseInt(line.get(1));
                        int OPLStype   = Integer.parseInt(line.get(2));
                        if ( OPLSatomTypes.size() != atomNumber - 1 )
                            throw new IllegalArgumentException("OPLS atom type out of order\n" + line.toString());
                        if ( OPLStype < 0 )
                            throw new IllegalArgumentException("can't have a negative OPLS atom type");
                        OPLSatomTypes.add(OPLStype);
                    }

                // read the surface tensions
                else if ( firstField.equals("surfacetension") )
                    {
                        if ( line.size() < 3 )
                            throw new IllegalArgumentException("unexpected number of fields for surface tension\n" + line.toString());
                        int atomNumber = Integer.parseInt(line.get(1));
                        double surfaceTension = Double.parseDouble(line.get(2));
                        if ( surfaceTensions.size() != atomNumber - 1 )
                            throw new IllegalArgumentException("surface tension out of order\n" + line.toString());
                        surfaceTensions.add(surfaceTension);
                    }
            }

        // check that we have all the information we need
        if ( OPLSatomTypes.size() != tempMolecule.contents.size() )
            throw new IllegalArgumentException(String.format("unexpected number of OPLS atom types (expected %d found %d)", tempMolecule.contents.size(), OPLSatomTypes.size()));
        if ( surfaceTensions.size() != tempMolecule.contents.size() )
            throw new IllegalArgumentException(String.format("unexpected number of surface tensions (expected %d found %d)", tempMolecule.contents.size(), surfaceTensions.size()));

        // create residue
        tempResidue = new Residue(aminoAcid, omega, phi, psi, chis,
                                  HN, N, O, C, CA, HA, description,
                                  prochiralConnection, atoms, isHairpin);

        // set atom types and surface tensions
        Map<Atom,Atom> atomMap = new HashMap<>();
        for (int i=0; i < atoms.size(); i++)
            {
                Atom a = atoms.get(i);
                int type2 = OPLSatomTypes.get(i);
                double surfaceTension = surfaceTensions.get(i);
                Atom newAtom = new Atom(a.element, a.position, a.type1, type2, surfaceTension);
                atomMap.put(a, newAtom);
            }
        tempResidue = tempResidue.moveAtoms(atomMap);

        tempNStickyConnection = movePair(tempNStickyConnection, atomMap);
        tempCStickyConnection = movePair(tempCStickyConnection, atomMap);
        tempMolecule          = tempMolecule.moveAtoms(atomMap);

        // set fields
        NStickyConnection = tempNStickyConnection;
        CStickyConnection = tempCStickyConnection;
        molecule          = tempMolecule;
        residue           = tempResidue;
    }
    
    /**
     * Parses lines of the form:<p>
     * field_title 1 2 3 4 [F]<p>
     * where 1 2 3 4 are Tinker atom types.
     * Will throw an error if there are multiple atoms corresponding to the same atom type.
     * @param molecule the molecule to grab the atoms from
     * @param line the newline-separated fields of the current line
     * @return the ProtoTorsion that corresponds to these atom types
     */
    private static ProtoTorsion parseProtoTorsion(Molecule molecule, List<String> line)
    {
        if ( line.get(1).toLowerCase().equals("null") )
            return null;
        int atom1type = Integer.parseInt(line.get(1));
        int atom2type = Integer.parseInt(line.get(2));
        int atom3type = Integer.parseInt(line.get(3));
        int atom4type = Integer.parseInt(line.get(4));
        Atom atom1    = getSingleAtom(molecule, atom1type);
        Atom atom2    = getSingleAtom(molecule, atom2type);
        Atom atom3    = getSingleAtom(molecule, atom3type);
        Atom atom4    = getSingleAtom(molecule, atom4type);
        return new ProtoTorsion(atom1, atom2, atom3, atom4);
    }

    /**
     * Parses lines of the form:<p>
     * field_title 1 2
     * where 1 2 are Tinker atom types.
     * Will throw an error if there are multiple atoms corresponding to the same atom type.
     * @param molecule the molecule to grab the atoms from
     * @param line the newline-separated fields of the current line
     * @return the pair that corresponds to these atom types
     */
    private static Pair<Atom,Atom> parsePair(Molecule molecule, List<String> line)
    {
        if ( line.get(1).toLowerCase().equals("null") )
            return null;
        int atom1type = Integer.parseInt(line.get(1));
        int atom2type = Integer.parseInt(line.get(2));
        Atom atom1    = getSingleAtom(molecule, atom1type);
        Atom atom2    = getSingleAtom(molecule, atom2type);
        return new Pair<Atom,Atom>(atom1, atom2);
    }

    /**
     * Parses lines of the form:<p>
     * field_title 1 2
     * where 1 2 are atom numbers 1,2,...,N
     * Will throw an error if the atom is not found.
     * @param molecule the molecule to grab the atoms from
     * @param line the newline-separated fields of the current line
     * @return the pair that corresponds to these atom types
     */
    private static Pair<Atom,Atom> parseAtomPair(Molecule molecule, List<String> line)
    {
        int atom1number = Integer.parseInt(line.get(1));
        int atom2number = Integer.parseInt(line.get(2));
        Atom atom1      = molecule.contents.get(atom1number-1);
        Atom atom2      = molecule.contents.get(atom2number-1);
        return new Pair<Atom,Atom>(atom1, atom2);
    }

    /**
     * Parses lines of the form:<p>
     * field_title 1 2
     * where 1 2 are atom numbers 1,2,...,N
     * Will throw an error if the atom is not found.
     * @param molecule the molecule to grab the atoms from
     * @param line the newline-separated fields of the current line
     * @return the pair that corresponds to these atom types
     */
    private static ProtoTorsion parseAtomQuadruple(Molecule molecule, List<String> line)
    {
        int atom1number = Integer.parseInt(line.get(1));
        int atom2number = Integer.parseInt(line.get(2));
        int atom3number = Integer.parseInt(line.get(3));
        int atom4number = Integer.parseInt(line.get(4));
        Atom atom1      = molecule.contents.get(atom1number-1);
        Atom atom2      = molecule.contents.get(atom2number-1);
        Atom atom3      = molecule.contents.get(atom3number-1);
        Atom atom4      = molecule.contents.get(atom4number-1);
        return new ProtoTorsion(atom1, atom2, atom3, atom4);
    }

    /**
     * Parses lines of the form:<p>
     * field_title 1 2 ...
     * where 1 2 ... are Tinker atom types.
     * @param molecule the molecule to grab the atoms from
     * @param line the newline-separated fields of the current line
     * @return the Atoms that corresponds to these atom types
     */
    private static List<Atom> parseMultipleAtoms(Molecule molecule, List<String> line)
    {
        if ( line.get(1).toLowerCase().equals("null") )
            return ImmutableList.of();
        List<Atom> returnList = new LinkedList<Atom>();
        for (int i=1; i < line.size(); i++)
            {
                String currentField = line.get(i);
                returnList.addAll(getMultipleAtoms(molecule, Integer.parseInt(currentField)));
            }
        return ImmutableList.copyOf(returnList);
    }

    /**
     * Parses lines of the form:<p>
     * field_title 1
     * where 1 is a Tinker atom types.
     * Will throw an error if there are multiple hits.
     * @param molecule the molecule to grab the atoms from
     * @param line the newline-separated fields of the current line
     * @param targetTypes the atom types we are looking for
     * @return the Atom that corresponds to this atom types
     */
    private static Atom parseSingleAtom(Molecule molecule, List<String> line, List<Integer> targetTypes)
    {
        if ( line.get(1).toLowerCase().equals("null") )
            return null;
        else if ( line.get(1).toLowerCase().equals("auto") )
            return getSingleAtom(molecule, targetTypes);
        else if ( line.get(2).toLowerCase().equals("a") )
            return molecule.contents.get( Integer.parseInt(line.get(1))-1 );
        return getSingleAtom(molecule, ImmutableList.of( Integer.parseInt(line.get(1)) ) );
    }

    // searches through the molecule for atoms matching this description
    private static List<Atom> getMultipleAtoms(Molecule molecule, int atomType)
    {
        List<Atom> results = new LinkedList<>();
        for ( Atom a : molecule.contents )
            {
                if ( a.type1 == atomType )
                    results.add(a);
            }
        return results;
    }

    // same as getAtom(int atomType), but only allows one result
    private static Atom getSingleAtom(Molecule molecule, List<Integer> targetTypes)
    {
        List<Atom> results = new ArrayList<>();
        for (Atom a : molecule.contents)
            if ( targetTypes.contains(a.type1) )
                results.add(a);
        if ( results.size() == 0 )
            throw new IllegalArgumentException("no matches found for atom types " + targetTypes.toString());
        else if ( results.size() > 1 )
            throw new IllegalArgumentException("multiple matches found for atom type " + targetTypes.toString());
        return results.get(0);
    }

    private static Atom getSingleAtom(Molecule molecule, int targetType)
    {
        return getSingleAtom(molecule, ImmutableList.of(targetType));
    }

    /**
     * Converts an atom pair using an atom map.
     */
    public static Pair<Atom,Atom> movePair(Pair<Atom,Atom> oldPair, Map<Atom,Atom> atomMap)
    {
        Atom tempAtom1 = oldPair.getFirst();
        Atom tempAtom2 = oldPair.getSecond();

        if ( atomMap.containsKey(tempAtom1) )
            tempAtom1 = atomMap.get(tempAtom1);
        else
            throw new IllegalArgumentException("first atom not found");

        if ( atomMap.containsKey(tempAtom2) )
            tempAtom2 = atomMap.get(tempAtom2);
        else
            throw new IllegalArgumentException("second atom not found");

        return new Pair<Atom,Atom>(tempAtom1, tempAtom2);
    }   

    /** for testing */
    public static void main(String[] args)
    {
        ProtoAminoAcidFile m = new ProtoAminoAcidFile("amino_acids/Ser.txt");
        ProtoAminoAcid protoAlanine = new ProtoAminoAcid(m);
        for (Atom a : m.molecule.contents)
            System.out.printf("%3d %s\n", m.molecule.contents.indexOf(a)+1, a.toFullString());
        System.out.println(protoAlanine);

    }
}
