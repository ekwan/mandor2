import java.util.*;
import java.io.*;
import org.apache.commons.math3.geometry.euclidean.threed.*;
import org.jgrapht.*;
import org.jgrapht.graph.*;
import org.jgrapht.alg.*;
import com.google.common.collect.*;
import java.util.concurrent.atomic.*;

/**
 * Represents an interaction in a molecule.
 */
public class Interaction implements Immutable, Comparable<Interaction>
{
    /** the atoms involved in this interaction */
    public final Set<Atom> atoms;

    /** the interaction energy */
    public final double interactionEnergy;

    /** the reference energies */
    public static final Map<AminoAcid,Double> REFERENCE_ENERGIES;

    /** static initializer */
    static
    {
        // read the OPLS reference energies from a file
        Map<AminoAcid,Double> tempMap = new TreeMap<AminoAcid,Double>();
        for (AminoAcid a : AminoAcid.values())
            tempMap.put(a, 0.0);
/*        OutputFileFormat file = new OutputFileFormat("amino_acids/oplsaa/reference.txt") {};
        for (List<String> line : file.fileContents)
            {
                String aaString = line.get(0);
                AminoAcid aminoAcid = AminoAcid.getAminoAcid(aaString);
                Double refEnergy = Double.valueOf(line.get(5));
                //Double refEnergy = 0.0;
                tempMap.put(aminoAcid, refEnergy);
            }*/
        REFERENCE_ENERGIES = ImmutableMap.copyOf(tempMap);
    }

    //public String description;

    /**
     * Constructor.
     * @param atoms the atoms involved in this interaction
     * @param interactionEnergy the interaction energy (kcal)
     */
    //public Interaction(Set<Atom> atoms, double interactionEnergy, String description)
    public Interaction(Set<Atom> atoms, double interactionEnergy)
    {
        this.atoms = ImmutableSet.copyOf(atoms);
        this.interactionEnergy = interactionEnergy;
        //this.description = description;
    }

    @Override
    public String toString()
    {
        return String.format("%s = %.4f kcal", atoms.toString(), interactionEnergy);
    }

    /**
     * Gives a nicer output with atom numbers as well.
     * @param molecule the molecule containing these atoms
     */
    public String toString(Molecule molecule)
    {
        List<String> description = new LinkedList<>();
        for (Atom a : atoms)
            description.add(molecule.getAtomString(a));
        return String.format("%s = %.4f kcal", description.toString(), interactionEnergy);
    }

    @Override
    public int hashCode()
    {
        return Objects.hash(atoms, interactionEnergy);
    }

    @Override
    public boolean equals(Object obj)
    {
        if ( obj == null )
            return false;
        if ( obj == this )
            return true;
        if ( !(obj instanceof Interaction) )
            return false;

        Interaction i = (Interaction)obj;
        if ( Objects.equals(atoms, i.atoms) &&
             interactionEnergy == i.interactionEnergy )
            return true;
        return false;
    }

    @Override
    public int compareTo(Interaction i)
    {
        return ComparisonChain.start().compare(i.interactionEnergy, interactionEnergy).result();
    }

    /** 
     * Utility method for getting the sidechain atoms of a peptide residue.
     * Could run into a problem if the peptide connectivity is concurrently modified!
     * @param peptide the peptide
     * @param residue a residue in the peptide
     * @return the sidechain atoms in the specified residue
     */
    public static Set<Atom> getSidechainAtoms(Peptide peptide, Residue residue)
    {
        if ( residue.isHairpin )
            return new HashSet<Atom>();
        
        Pair<Atom,Atom> prochiralConnection = residue.prochiralConnection;
        Set<Atom> sidechainAtoms = null;
        if (residue.aminoAcid.isProline())
            {
                // temporarily disconnect atom2 and atom3 of chi4
                ProtoTorsion chi3 = residue.chis.get(2);
                Atom atom3 = chi3.atom3;
                Atom atom4 = chi3.atom4;
                SimpleWeightedGraph<Atom, DefaultWeightedEdge> connectivity = peptide.connectivity;
                connectivity.removeEdge(atom3, atom4);
                sidechainAtoms = peptide.getHalfGraph(prochiralConnection.getFirst(),prochiralConnection.getSecond());
                DefaultWeightedEdge e = connectivity.addEdge(atom3, atom4);
                connectivity.setEdgeWeight(e, 1.0);
            }
        else
            {
                sidechainAtoms = peptide.getHalfGraph(prochiralConnection.getFirst(),prochiralConnection.getSecond());
                sidechainAtoms = new HashSet<Atom>(sidechainAtoms);
                sidechainAtoms.add(residue.backboneHN);
            }
        return sidechainAtoms;
    }

    //public static AtomicInteger counter = new AtomicInteger();

    /**
     * Creates a matrix of rotamer energies for a particular peptide.  Diagonal entries are rotamer self-energies,
     * which are defined as the interactions inside the rotamer as well as rotamer-backbone interactions.  The
     * entries are arranged in sequence order (0, 1, ..., n-1).  For example, entry (0,1) (or (1,0), since the matrix
     * is symmetric) corresponds to the interaction energy between the rotamer at position 0 and the rotamer at position 1.
     * The backbone is treated as the n-th residue, such that the (n,n) entry is the backbone self-energy.
     * @param peptide the peptide that contains the rotamers
     * @param interactions the interactions in this peptide
     */
    public static Double[][] getRotamerEnergyMatrix(Peptide peptide, List<Interaction> interactions)
    {
        // initialize matrix where the results will go
        // the backbone is treated as the n+1-th residue,
        // where n is the number of residues
        int numberOfResidues = peptide.sequence.size();
        double[][] energyMatrix = new double[numberOfResidues+1][numberOfResidues+1];

        // get all rotamer atoms
        List<Set<Atom>> rotamerAtoms = new ArrayList<>(numberOfResidues);
        Set<Atom> allRotamerAtoms = new HashSet<>();
        for (Residue r : peptide.sequence)
            {
                Set<Atom> atoms = getSidechainAtoms(peptide,r);
                rotamerAtoms.add(atoms);
                allRotamerAtoms.addAll(atoms);
            }

        // get backbone atoms
        Set<Atom> backboneAtoms = new HashSet<>();
        for (Atom a : peptide.contents)
            {
                if ( allRotamerAtoms.contains(a) )
                    continue;
                backboneAtoms.add(a);
            }

        // classify the interactions
        /*
        int backboneInteractions = 0;
        double backboneEnergy = 0.0;

        int count = counter.getAndIncrement();
        
        String debugFilename = String.format("test_peptides/peptide_%05d.gjf", count);
        GaussianInputFile f = new GaussianInputFile(peptide);
        f.write(debugFilename);

        String debugFilename2 = String.format("test_peptides/backbone_%05d.gjf", count);
        f = new GaussianInputFile(peptide, backboneAtoms);
        f.write(debugFilename2);

        ArrayList<Atom> backboneList = new ArrayList<Atom>(backboneAtoms);
        Collections.sort(backboneList);

        List<String> backboneDescriptionList = new ArrayList<>();
        */
        for (Interaction interaction : interactions)
            {
                List<Atom> backboneAtomsList = ImmutableList.copyOf(interaction.atoms);

                // figure out which groups this interaction belongs to
                boolean inBackbone = false;
                boolean inRotamer1 = false;
                int rotamer1index  = -1;
                boolean inRotamer2 = false;
                int rotamer2index  = -1;

                for (Atom a : interaction.atoms)
                    {
                        if ( backboneAtoms.contains(a) )
                            inBackbone = true;
                        else
                            {
                                for (int rotamerIndex = 0; rotamerIndex < rotamerAtoms.size(); rotamerIndex++)
                                    {
                                        Set<Atom> thisRotamerAtoms = rotamerAtoms.get(rotamerIndex);
                                        if ( thisRotamerAtoms.contains(a) )
                                            {
                                                if ( rotamerIndex == rotamer1index || rotamerIndex == rotamer2index )
                                                    break;
                                                else if ( rotamer1index == -1 )
                                                    {
                                                        inRotamer1 = true;
                                                        rotamer1index = rotamerIndex;
                                                        break;
                                                    }
                                                else if ( rotamer2index == -1 )
                                                    {
                                                        inRotamer2 = true;
                                                        rotamer2index = rotamerIndex;
                                                        break;
                                                    }
                                                else
                                                    throw new IllegalArgumentException("interaction cannot be in three rotamers");
                                            }
                                    }
                            }
                    }

                // put the energy in the correct bin
                double interactionEnergy = interaction.interactionEnergy;
                if      (  inBackbone && !inRotamer1 && !inRotamer2 )
                    {
                        //System.out.println("backbone only");
                        energyMatrix[numberOfResidues][numberOfResidues] += interactionEnergy;
                        /*
                        backboneInteractions++;
                        backboneEnergy += interactionEnergy;
                        
                        int number1 = backboneList.indexOf(backboneAtomsList.get(0))+1;
                        int number2 = backboneList.indexOf(backboneAtomsList.get(1))+1;
                        if ( number1 < number2 )
                            backboneDescriptionList.add(String.format("\n%d-%d  :  %s", number1, number2, interaction.description));
                        else
                            backboneDescriptionList.add(String.format("\n%d-%d  :  %s", number2, number1, interaction.description));
                        */
                    }
                else if (  inBackbone &&  inRotamer1 && !inRotamer2 )
                    {
                        //System.out.println("position " + rotamer1index + " backbone");
                        energyMatrix[rotamer1index][numberOfResidues] += interactionEnergy;
                        energyMatrix[numberOfResidues][rotamer1index] += interactionEnergy;
                    }
                else if (  inBackbone && !inRotamer1 &&  inRotamer2 )
                    {
                        //System.out.println("position " + rotamer2index + " backbone");
                        energyMatrix[rotamer2index][numberOfResidues] += interactionEnergy;
                        energyMatrix[numberOfResidues][rotamer2index] += interactionEnergy;
                    }
                else if ( !inBackbone &&  inRotamer1 && !inRotamer2 )
                    {
                        //System.out.println("position " + rotamer1index + " only");
                        energyMatrix[rotamer1index][rotamer1index] += interactionEnergy;
                    }
                else if ( !inBackbone && !inRotamer1 &&  inRotamer2 )
                    {
                        //System.out.println("position " + rotamer2index + " only");
                        energyMatrix[rotamer2index][rotamer2index] += interactionEnergy;
                    }
                else if ( !inBackbone &&  inRotamer1 &&  inRotamer2 )
                    {
                        //System.out.println("position " + rotamer1index + ", " + rotamer2index + " interaction");
                        energyMatrix[rotamer1index][rotamer2index] += interactionEnergy;
                        energyMatrix[rotamer2index][rotamer1index] += interactionEnergy;
                    }
                else
                    throw new IllegalArgumentException("error assigning interaction");
            }
        
        /*System.out.printf("[%3d] %d backbone self-interactions, %d backbone atoms, %.4f kcal\n", count, backboneInteractions, backboneAtoms.size(), backboneEnergy);

        String debugFilename3 = String.format("test_peptides/backbone_%05d.txt", count);
        
        String backboneDescription = "";
        Collections.sort(backboneDescriptionList);
        for (String s : backboneDescriptionList)
            backboneDescription += s;
        InputFileFormat.writeStringToDisk(backboneDescription, debugFilename3);
        */

        // return result
        Double[][] resultMatrix = new Double[numberOfResidues+1][numberOfResidues+1];
        for (int i=0; i < numberOfResidues+1; i++)
            {
                for (int j=0; j < numberOfResidues+1; j++)
                    {
                        resultMatrix[i][j] = Double.valueOf(energyMatrix[i][j]);
                    }
            }

        return resultMatrix;
    }

    public static double getReferenceEnergy(Rotamer rotamer)
    {
        AminoAcid aminoAcid = rotamer.aminoAcid;
        Double referenceEnergy = REFERENCE_ENERGIES.get(aminoAcid);
        if ( referenceEnergy == null )
            throw new NullPointerException("null reference energy for " + aminoAcid.toString());
        return referenceEnergy;
    }

    /**
     * Returns the energy of the specified rotamer.  This includes the energy of the rotamer itself and its interactions
     * with the backbone.  Includes reference energies if requested.
     * @param rotamer the rotamer to get the energy of
     * @param energyMatrix see {@link #getRotamerEnergyMatrix}
     * @param useReferenceEnergies if true use the OPLS reference energies
     * @return the self-energy of this rotamer
     */
    public static Double getRotamerEnergy(Rotamer rotamer, Double[][] energyMatrix, boolean useReferenceEnergies )
    {
        double rawEnergy = energyMatrix[rotamer.sequenceIndex][energyMatrix.length-1] + energyMatrix[rotamer.sequenceIndex][rotamer.sequenceIndex];
        if ( ! useReferenceEnergies )
            return rawEnergy;
        return rawEnergy - getReferenceEnergy(rotamer);
    }

    /**
     * Get rotamer energies with reference energies included.
     * See {@link #getRotamerEnergy(Rotamer,Double[][],boolean)}.
     * @return the energy of the rotamer
     */
    public static Double getRotamerEnergy(Rotamer rotamer, Double[][] energyMatrix)
    {
        return getRotamerEnergy(rotamer, energyMatrix, true);
    }

    /**
     * Returns the energy of the specified rotamer pair.
     * @param rotamer1 a rotamer in the pair
     * @param rotamer2 the other rotamer in the pair
     * @param energyMatrix see {@link #getRotamerEnergyMatrix}
     */
    public static Double getRotamerInteractionEnergy(Rotamer rotamer1, Rotamer rotamer2, Double[][] energyMatrix)
    {
        int row = rotamer1.sequenceIndex;
        int col = rotamer2.sequenceIndex;
        return energyMatrix[row][col];
    }

    /**
     * Returns the energy of the backbone.
     * @param energyMatrix see {@link #getRotamerEnergyMatrix}
     */
    public static Double getBackboneEnergy(Double[][] energyMatrix)
    {
        return energyMatrix[energyMatrix.length-1][energyMatrix.length-1];
    }

    /** for testing */
    public static void main(String[] args)
    {
    }
}
