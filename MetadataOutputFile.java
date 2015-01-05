import java.io.*;
import java.util.*;
import com.google.common.collect.*;

/**
 * Metadata files contain information about special atoms, torsions,
 * , etc for ProtoAminoAcids.  toString, equals, and hashCode not implemented.
 * This class is immutable.
 */
public class MetadataOutputFile extends OutputFileFormat implements Serializable, Immutable
{
    /** for serialization */
    public static final long serialVersionUID = 1L;

    // see Residue.java documentation
    public final AminoAcid aminoAcid;
    public final ResidueType residueType;
    public final ProtoTorsion omega;
    public final ProtoTorsion phi;
    public final ProtoTorsion psi;
    public final List<ProtoTorsion> chis;
    public final Atom backboneHN;
    public final Atom imidazoleHN;
    public final Atom imidazoleN;
    public final List<Atom> otherHNatoms;
    public final Pair<Atom,Atom> NStickyConnection;
    public final Pair<Atom,Atom> CStickyConnection;
    public final Pair<Atom,Atom> prochiralConnection;
    public final double referenceEnergy;
    public final String description;

    /** the geometry in the corresponding Tinker XYZ file */
    public final Molecule molecule;

    /** the Tinker XYZ file that contains the gometry for this amino acid */
    public final TinkerXYZOutputFile geometryFile;

    /** filename of the Tinker XYZ file that contains the geometry for this amino acid */
    public final String XYZfilename;

    /** public constructor */
    public MetadataOutputFile(String filename)
    {
        super(filename);

        // check that this is a valid file
        String headerString = stringRepresentation.split("\n")[0];
        if ( headerString.indexOf("Amino Acid Metadata") == -1 )
            throw new IllegalArgumentException("this does not appear to be an amino acid metadata file");

        // create temporary variables
        // required, otherwise, we'd have to try and set the final variables inside the loop
        // ...which isn't allowed
        AminoAcid tempAminoAcid = null;
        ResidueType tempResidueType = null;
        ProtoTorsion tempOmega = null;
        ProtoTorsion tempPhi = null;
        ProtoTorsion tempPsi = null;
        List<ProtoTorsion> tempChis = new LinkedList<>();
        Atom tempBackboneHN = null;
        Atom tempImidazoleHN = null;
        Atom tempImidazoleN = null;
        List<Atom> tempOtherHNatoms = new LinkedList<>();
        Pair<Atom,Atom> tempNStickyConnection = null;
        Pair<Atom,Atom> tempCStickyConnection = null;
        Pair<Atom,Atom> tempProchiralConnection = null;
        Double tempReferenceEnergy = null;
        String tempDescription = "";

        Molecule tempMolecule = null;
        TinkerXYZOutputFile tempGeometryFile = null;
        String tempXYZfilename = null;

        // parse the things that don't depend on the molecular specification
        for (List<String> line : fileContents)
            {
                // check for blank or comment lines
                if ( line.size() < 2 || line.get(0).startsWith("!") )
                    continue;

                String firstField = line.get(0).toLowerCase();
                
                if ( firstField.equals("amino_acid") )
                    tempAminoAcid = AminoAcid.getAminoAcid(line.get(1));
                else if ( firstField.equals("xyz_file") )
                    tempXYZfilename = line.get(1);
                else if ( firstField.equals("residue_type") )
                    tempResidueType = ResidueType.getResidueType(line.get(1));
                else if ( firstField.equals("reference_energy") )
                    tempReferenceEnergy = Double.parseDouble(line.get(1));
                else if ( firstField.equals("description") )
                    {
                        for ( int i=1; i < line.size(); i++ )
                            {
                                tempDescription = tempDescription + line.get(i);
                                if ( i < line.size() )
                                    tempDescription = tempDescription + " ";
                            }
                    }
            }

        // read the molecular specification
        XYZfilename = tempXYZfilename;
        geometryFile = new TinkerXYZOutputFile(Settings.AMINO_ACID_DIRECTORY + XYZfilename);
        molecule = geometryFile.molecule;

        // set all geometry-related fields
        for (List<String> line : fileContents)
            {
                // check for blank or comment lines
                if ( line.size() < 2 || line.get(0).startsWith("!") )
                    continue;
                String firstField = line.get(0).toLowerCase();
               
                // parse fields by reading the Tinker atom type
                // and locating all the atoms that match that description
                // in the current molecule
                if ( firstField.equals("phi") )
                    {
                        tempPhi = parseProtoTorsion(line);
                        if ( tempPhi == null )
                            throw new NullPointerException("phi cannot be null");
                        if ( line.size() >= 6 && line.get(5).toUpperCase().equals("F") )
                            tempFrozenTorsions.add(tempPhi);
                    }
                else if ( firstField.equals("psi") )
                    {
                        tempPsi = parseProtoTorsion(line);
                        if ( tempPsi == null )
                            throw new NullPointerException("psi cannot be null");
                        if ( line.size() >= 6 && line.get(5).toUpperCase().equals("F") )
                            tempFrozenTorsions.add(tempPsi);
                    }
                else if ( firstField.equals("omega") )
                    {
                        tempOmega = parseProtoTorsion(line);
                        if ( tempOmega == null )
                            throw new NullPointerException("omega cannot be null");
                        if ( line.size() >= 6 && line.get(5).toUpperCase().equals("F") )
                            tempFrozenTorsions.add(tempOmega);
                    }
               else if ( firstField.startsWith("chi") )
                    {
                        ProtoTorsion thisTorsion = null;
                        if ( line.size() >= 6 && ( line.get(5).toUpperCase().equals("A") || line.get(5).toUpperCase().equals("AF") ) )
                            thisTorsion = parseAtomQuadruple(line);
                        else
                            thisTorsion = parseProtoTorsion(line);

                        if ( thisTorsion != null )
                            {
                                tempChis.add(thisTorsion);
                                int expectedSize = Integer.parseInt(firstField.replaceAll("[^0-9]",""));
                                if ( tempChis.size() != expectedSize )
                                    throw new IllegalArgumentException("are the chis out of order?");
                                if ( line.size() >= 6 && ( line.get(5).toUpperCase().equals("F") || line.get(5).toUpperCase().equals("AF") ) )
                                    tempFrozenTorsions.add(thisTorsion);
                            }
                    }
                else if ( firstField.equals("backbone_hn") )
                    tempBackboneHN = parseSingleAtom(line);
                else if ( firstField.equals("imidazole_hn") )
                    tempImidazoleHN = parseSingleAtom(line);
                else if ( firstField.equals("imidazole_n") )
                    tempImidazoleN = parseSingleAtom(line);
                else if ( firstField.equals("other_hn") )
                    tempOtherHNatoms = parseMultipleAtoms(line);
                else if ( firstField.equals("nstickyconnection") )
                    tempNStickyConnection = parsePair(line);
                else if ( firstField.equals("cstickyconnection") )
                    tempCStickyConnection = parsePair(line);
                else if ( firstField.equals("prochiral_connection") )
                    {
                        if ( line.size() >= 4 && line.get(3).toUpperCase().equals("A") )
                            tempProchiralConnection = parseAtomPair(line);
                        else
                            tempProchiralConnection = parsePair(line);
                    }
            }
        
        // consistency check
        if ( tempNStickyConnection == null )
            throw new NullPointerException("N-sticky connection is required!");
        if ( tempNStickyConnection.getFirst().element != Element.NITROGEN ||
             tempNStickyConnection.getSecond().element != Element.CARBON )
            throw new IllegalArgumentException("Illegal atom type for NStickyConnection");
        
        if ( tempCStickyConnection == null )
            throw new NullPointerException("N-sticky connection is required!");
        if ( tempCStickyConnection.getFirst().element != Element.CARBON ||
             tempCStickyConnection.getSecond().element != Element.NITROGEN )
            throw new IllegalArgumentException("Illegal atom type for CStickyConnection");

        if ( tempAminoAcid == AminoAcid.HIS )
            {
                if ( tempImidazoleHN == null || tempImidazoleN == null )
                    throw new IllegalArgumentException("histidine HN and N must be specified!");
            }
        else
            {
                if ( tempImidazoleHN != null || tempImidazoleN != null )
                    throw new IllegalArgumentException("should not specify imidazole atoms for a non-histidine");
            }
        
        if ( tempBackboneHN == null && !tempAminoAcid.isProline() )
            throw new NullPointerException("backbone HN must be specified");

        if ( tempReferenceEnergy == null )
            throw new NullPointerException("reference energy must be set");

        if ( tempProchiralConnection == null )
            throw new NullPointerException("must specify a prochiral connection!");

        // clean up
        aminoAcid = tempAminoAcid;
        residueType = tempResidueType;
        omega = tempOmega;
        phi = tempPhi;
        psi = tempPsi;
        chis = ImmutableList.copyOf(tempChis);
        backboneHN = tempBackboneHN;
        imidazoleHN = tempImidazoleHN;
        imidazoleN = tempImidazoleN;
        otherHNatoms = ImmutableList.copyOf(tempOtherHNatoms);
        NStickyConnection = tempNStickyConnection;
        CStickyConnection = tempCStickyConnection;
        prochiralConnection = tempProchiralConnection;
        referenceEnergy = tempReferenceEnergy;
        description = tempDescription;
    }

    /**
     * Parses lines of the form:<p>
     * field_title 1 2 3 4 [F]<p>
     * where 1 2 3 4 are Tinker atom types.
     * Will throw an error if there are multiple atoms corresponding to the same atom type.
     * @param line the newline-separated fields of the current line
     * @return the ProtoTorsion that corresponds to these atom types
     */
    private ProtoTorsion parseProtoTorsion(List<String> line)
    {
        if ( line.get(1).toLowerCase().equals("null") )
            return null;
        int atom1type = Integer.parseInt(line.get(1));
        int atom2type = Integer.parseInt(line.get(2));
        int atom3type = Integer.parseInt(line.get(3));
        int atom4type = Integer.parseInt(line.get(4));
        Atom atom1    = getSingleAtom(atom1type);
        Atom atom2    = getSingleAtom(atom2type);
        Atom atom3    = getSingleAtom(atom3type);
        Atom atom4    = getSingleAtom(atom4type);
        return new ProtoTorsion(atom1, atom2, atom3, atom4);
    }

    /**
     * Parses lines of the form:<p>
     * field_title 1 2
     * where 1 2 are Tinker atom types.
     * Will throw an error if there are multiple atoms corresponding to the same atom type.
     * @param line the newline-separated fields of the current line
     * @return the pair that corresponds to these atom types
     */
    private Pair<Atom,Atom> parsePair(List<String> line)
    {
        if ( line.get(1).toLowerCase().equals("null") )
            return null;
        int atom1type = Integer.parseInt(line.get(1));
        int atom2type = Integer.parseInt(line.get(2));
        Atom atom1    = getSingleAtom(atom1type);
        Atom atom2    = getSingleAtom(atom2type);
        return new Pair<Atom,Atom>(atom1, atom2);
    }

    /**
     * Parses lines of the form:<p>
     * field_title 1 2
     * where 1 2 are atom numbers 1,2,...,N
     * Will throw an error if the atom is not found.
     * @param line the newline-separated fields of the current line
     * @return the pair that corresponds to these atom types
     */
    private Pair<Atom,Atom> parseAtomPair(List<String> line)
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
     * @param line the newline-separated fields of the current line
     * @return the pair that corresponds to these atom types
     */
    private ProtoTorsion parseAtomQuadruple(List<String> line)
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
     * @param line the newline-separated fields of the current line
     * @return the Atoms that corresponds to these atom types
     */
    private List<Atom> parseMultipleAtoms(List<String> line)
    {
        if ( line.get(1).toLowerCase().equals("null") )
            return ImmutableList.of();
        List<Atom> returnList = new LinkedList<Atom>();
        for (int i=1; i < line.size(); i++)
            {
                String currentField = line.get(i);
                if ( currentField.startsWith("!") )
                    break;
                returnList.addAll(getMultipleAtoms(Integer.parseInt(currentField)));
            }
        return ImmutableList.copyOf(returnList);
    }

    /**
     * Parses lines of the form:<p>
     * field_title 1
     * where 1 is a Tinker atom types.
     * Will throw an error if there are multiple hits.
     * @param line the newline-separated fields of the current line
     * @return the Atom that corresponds to this atom types
     */
    private Atom parseSingleAtom(List<String> line)
    {
        if ( line.get(1).toLowerCase().equals("null") )
            return null;
        return getSingleAtom(Integer.parseInt(line.get(1)));
    }

    // searches through the molecule for atoms matching this description
    private List<Atom> getMultipleAtoms(int atomType)
    {
        List<Atom> results = new LinkedList<>();
        for ( Atom a : molecule.contents )
            {
                if ( a.tinkerAtomType == atomType )
                    results.add(a);
            }
        return results;
    }

    // same as getAtom(int atomType), but only allows one result
    private Atom getSingleAtom(int atomType)
    {
        List<Atom> results = getMultipleAtoms(atomType);
        if ( results.size() == 0 )
            throw new IllegalArgumentException("no matches found for atom type " + atomType);
        else if ( results.size() > 1 )
            throw new IllegalArgumentException("multiple matches found for atom type " + atomType);
        return results.get(0);
    }

    /** for testing */
    public static void main(String[] args)
    {
        MetadataOutputFile m = new MetadataOutputFile("amino_acids/Tyr.txt");
        ProtoAminoAcid protoTyrosine = new ProtoAminoAcid(m);
        System.out.println(protoTyrosine);
        InputFileFormat.writeStringToDisk(protoTyrosine.molecule.toGaussianString(), "amino_acids/Tyr.gjf");
    }
}
