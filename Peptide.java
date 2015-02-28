import java.util.*;
import java.io.*;
import com.google.common.collect.*;
import org.apache.commons.math3.geometry.euclidean.threed.*;
import org.jgrapht.*;
import org.jgrapht.graph.*;
import org.jgrapht.alg.*;
import java.util.concurrent.*;

/**
 * This class represents a peptide.  It is immutable.
 */
public class Peptide extends Molecule implements Immutable, Serializable, Comparable<Peptide>
{
    /** For serialization. */
    public static final long serialVersionUID = 1L;
    
    /** The residues, from the N-terminus to the C-terminus. */
    public final List<Residue> sequence;
    
    /** The energy broken down by residue. */
    public final EnergyBreakdown energyBreakdown;

    /** Creates a Peptide. */
    public Peptide(String name, List<Atom> contents, SimpleWeightedGraph<Atom,DefaultWeightedEdge> connectivity,
                   List<Residue> sequence, EnergyBreakdown energyBreakdown)
    {
        super(name,contents,connectivity);
        this.sequence = ImmutableList.copyOf(sequence);
        this.energyBreakdown = energyBreakdown;
    }

    /**
     * Factory method to create a peptide given a map of old atoms to new atoms.  Should be used
     * to move atoms.
     * @param atomMap a map from old atoms to new atoms (does not have to include all atoms)
     * @return the new Peptide
     */
    public Peptide moveAtoms2(Map<Atom,Atom> atomMap)
    {
        Molecule oldMolecule = (Molecule)this;
        Molecule newMolecule = oldMolecule.moveAtoms(atomMap);

        // create new sequence/residues
        List<Residue> newSequence = new ArrayList<>(sequence.size());
        for (Residue r : sequence)
            newSequence.add(r.moveAtoms(atomMap));
        newSequence = ImmutableList.copyOf(newSequence);
        
        // return result
        // assumes the new peptide has the same name and won't have an EnergyBreakdown yet
        return new Peptide(newMolecule.name, newMolecule.contents, newMolecule.connectivity, newSequence, EnergyBreakdown.BLANK);
    }

    /**
     * Returns a new Peptide that has the same fields as this one, except for the name.
     * @param newName the new name
     * @return the new Peptide
     */
    public Peptide setName(String newName)
    {
        return new Peptide(newName, contents, connectivity, sequence, energyBreakdown);
    }

    /**
     * Returns a new Peptide that has the same fields as this one, but with a different EnergyBreakdown.
     * @param newEnergyBreakdown the new energy breakdown to use
     * @return the new Peptide
     */
    public Peptide setEnergyBreakdown(EnergyBreakdown newEnergyBreakdown)
    {
        //if ( energyBreakdown != null && energyBreakdown != EnergyBreakdown.BLANK )
        //    System.out.println("warning: set energy breakdown when not null");
        String newName = name.split("@")[0] + String.format(" @ %.2f kcal", newEnergyBreakdown.totalEnergy);
        return new Peptide(newName, contents, connectivity, sequence, newEnergyBreakdown);
    }

    /**
     * Changes the geometry of this Peptide.
     * @param newMolecule the new molecule geometry to use
     * @return the new Peptide
     */
    public Peptide setMolecule(Molecule newMolecule)
    {
        if ( contents.size() != newMolecule.contents.size() )
            throw new IllegalArgumentException("contents size mismatch");
        Molecule oldMolecule = (Molecule)this;
        Map<Atom,Atom> atomMap = getAtomMap(oldMolecule,newMolecule);
        return moveAtoms2(atomMap);
    }

    /**
     * Changes the geometry of this Peptide.
     * @param newMolecule the new molecule geometry to use
     * @return the new Peptide
     */
    public Peptide setPositions(Molecule newMolecule)
    {
        if ( contents.size() != newMolecule.contents.size() )
            throw new IllegalArgumentException("contents size mismatch");
        Map<Atom,Atom> atomMap = new HashMap<>();
        for (int i=0; i < contents.size(); i++)
            {
                Atom oldAtom = contents.get(i);
                Vector3D newPosition = newMolecule.contents.get(i).position;
                Atom newAtom = oldAtom.moveAtom(newPosition);
                atomMap.put(oldAtom, newAtom);
            }
        return moveAtoms2(atomMap);
    }

    /**
     * Changes the geometry of this Peptide.
     * @param positions the new geometry to use
     * @return the new Peptide
     */
    public Peptide setMolecule(List<Vector3D> positions)
    {
        if ( contents.size() != positions.size() )
            throw new IllegalArgumentException("size mismatch");
        Map<Atom,Atom> atomMap = new HashMap<>();
        for (int i=0; i < contents.size(); i++)
            {
                Atom oldAtom = contents.get(i);
                Atom newAtom = oldAtom.moveAtom(positions.get(i));
                atomMap.put(oldAtom,newAtom);
            }
        return moveAtoms2(atomMap);
    }

    /**
     * Provides an atom map from one molecule to another.  This assumes
     * molecule1 and molecule2 are different conformations.
     * @param molecule1 the from atoms
     * @param molecule2 the to atoms
     * @return a map from atoms in molecule1 to atoms in molecule2 (won't include atoms that don't move)
     */
    public static Map<Atom,Atom> getAtomMap(Molecule molecule1, Molecule molecule2)
    {
        if ( molecule1.contents.size() != molecule2.contents.size() )
            throw new IllegalArgumentException("atom list size mismatch");
        Map<Atom,Atom> returnMap = new HashMap<>();
        for (int i=0; i < molecule1.contents.size(); i++)
            {
                Atom a1 = molecule1.contents.get(i);
                Atom a2 = molecule2.contents.get(i);
                if ( ! a1.equals(a2) )
                    returnMap.put(a1,a2);
            }
        return returnMap;
    }

    /**
     * Writes the specified peptides to disk in gjf format.
     * @param peptides the peptides to write to disk
     * @param prefix the location and prefix to send the peptides to (e.g., "test_peptides/prefix_")
     * @param digits how many digits to use in the filename
     * @param maxNumber the maximum number of files to write
     */
    public static void writeGJFs(List<Peptide> peptides, String prefix, int digits, int maxNumber)
    {
        if ( peptides == null )
            return;
        if ( digits < 1 )
            throw new IllegalArgumentException("must use at least one digit");
        for (int i=0; i < Math.min(maxNumber,peptides.size()); i++)
            {
                Peptide peptide = peptides.get(i);
                String filename = String.format("%s%0" + digits + "d.gjf", prefix, i);
                GaussianInputFile f = new GaussianInputFile(peptide);
                f.write(filename);
            }
    }

    /**
     * Writes the specified peptides to disk in gjf format.
     * @param peptides the peptides to write to disk
     * @param prefix the location and prefix to send the peptides to (e.g., "test_peptides/prefix_")
     * @param digits how many digits to use in the filename
     * @param digits how many digits to use in the filename
     */
    public static void writeCHKs(List<Peptide> peptides, String prefix, int digits, int maxNumber)
    {
        if ( peptides == null )
            return;
        if ( digits < 1 )
            throw new IllegalArgumentException("must use at least one digit");
        for (int i=0; i < Math.min(maxNumber,peptides.size()); i++)
            {
                Peptide peptide = peptides.get(i);
                String filename = String.format("%s%0" + digits + "d.chk", prefix, i);
                peptide.checkpoint(filename);
            }
    }

    /**
     * Looks for pairs of backbone atoms that should be checked for clashes.
     * Only backbone atoms that are more than two bonds apart are checked.
     * Hairpin positions are treated as part of the backbone.  All other positions
     * are treated as having rotamers, even if it's just the H of glycine.
     * HNs are treated as part of the backbone for this method.
     * @return pairs of indices of backbone atoms that should be checked for clashes
     */
    public List<Pair<Integer,Integer>> getBackbonePairs()
    {
        // get all rotamer atoms
        int numberOfResidues = sequence.size();
        Set<Atom> allRotamerAtoms = new HashSet<>();
        for (Residue r : sequence)
            {
                if ( r.isHairpin )
                    continue;
                allRotamerAtoms.addAll( RotamerFactory.getSidechainAtoms(this,r,false) );
            }

        // get backbone atoms
        Set<Atom> backboneAtoms = new HashSet<>();
        for (Atom a : contents)
            {
                if ( allRotamerAtoms.contains(a) )
                    continue;
                backboneAtoms.add(a);
            }
        if ( backboneAtoms.size() == 0 )
            throw new IllegalArgumentException("expected to find some backbone atoms");
        
        // create all pairs
        List<Atom> backboneAtomsList = new ArrayList<>(backboneAtoms);
        int expectedSize = backboneAtomsList.size() * (backboneAtomsList.size() - 1);
        List<Pair<Integer,Integer>> returnList = new ArrayList<>(expectedSize);
        for (int i=0; i < backboneAtomsList.size()-1; i++)
            {
                Atom atomI = contents.get(i);
                for (int j=i+1; j < backboneAtomsList.size(); j++)
                    {
                        Atom atomJ = contents.get(j);
                        if ( !areSeparated(atomI,atomJ) )
                            continue;
                        Pair<Integer,Integer> pair = new Pair<>(i,j);
                        returnList.add(pair);
                    }
            }

        // return result
        return ImmutableList.copyOf(returnList);
    }

    /**
     * Checks if the backbone atoms in the specified peptide clash. 
     * @param backbonePairs pairs of indices of backbone atoms that should be checked for clashes
     * @return true if there is a clash
     */
    public boolean hasBackboneClash(List<Pair<Integer,Integer>> backbonePairs)
    {
        for (Pair<Integer,Integer> pair : backbonePairs)
            {
                int i = pair.getFirst();
                int j = pair.getSecond();
                Atom atomI = contents.get(i);
                Atom atomJ = contents.get(j);
                if ( RotamerSpace.clashes(atomI,atomJ) )
                    return true;
            }
        return false;
    }

    /**
     * Writes the peptide to disk.  Will not die from exceptions.
     * @param filename the filename to write to
     */
    public void checkpoint(String filename)
    {
        if ( filename == null || filename.length() == 0 )
            throw new NullPointerException("must specify a filename");
        try
            {
                 FileOutputStream fileOut = new FileOutputStream(filename);
                 ObjectOutputStream out = new ObjectOutputStream(fileOut);
                 out.writeObject(this);
                 out.close();
                 fileOut.close();
            }
        catch (Exception e)
            {
                e.printStackTrace();
            }
    }

    /**
     * Loads a peptide checkpoint from disk.  WIll not die from exceptions.
     * @param filename the filename to write to
     * @return the deserialized peptide
     */
    public static Peptide load(String filename)
    {
        if ( filename == null || filename.length() == 0 )
            throw new NullPointerException("must specify a filename");
        try
            {
                FileInputStream fileIn = new FileInputStream(filename);
                ObjectInputStream in = new ObjectInputStream(fileIn);
                Peptide p = (Peptide)in.readObject();
                in.close();
                fileIn.close();
                return p;
            }
        catch (Exception e)
            {
                e.printStackTrace();
                return null;
            }
    }

    /**
     * Compares peptides on the basis of their total energy in
     * energy breakdown (allows for sorting lists in ascending order) 
     * @param p2 the peptide to compare this peptide to
     */
    @Override
    public int compareTo(Peptide p2)
    {
        if ( this.energyBreakdown == null || p2.energyBreakdown == null ||
             this.energyBreakdown == EnergyBreakdown.BLANK || p2.energyBreakdown == EnergyBreakdown.BLANK )
            throw new NullPointerException("null or blank energy breakdown -- cannot compare");
        if ( this.energyBreakdown.type != p2.energyBreakdown.type )
            throw new IllegalArgumentException("energy type mismatch");
        return energyBreakdown.totalEnergy > p2.energyBreakdown.totalEnergy ? 1 : (energyBreakdown.totalEnergy < p2.energyBreakdown.totalEnergy ? -1 : 0);
    }

    @Override
    public int hashCode()
    {
        return Objects.hash(name, contents, connectivity, sequence, energyBreakdown);
    }

    @Override
    public boolean equals(Object obj)
    {
        if ( obj == null )
            return false;
        if ( obj == this )
            return true;
        if ( !(obj instanceof Peptide) )
            return false;

        Peptide p = (Peptide)obj;
        if ( name.equals(p.name) &&
             contents.equals(p.contents) &&
             connectivity.equals(p.connectivity) &&
             sequence.equals(p.sequence) &&
             Objects.equals(energyBreakdown, p.energyBreakdown) )
            return true;
        return false;
    }
}
