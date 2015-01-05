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
    public static final long serialVersionUID = 1L;

    /** the sequence of residues in a peptide starting from the N terminal to the C terminal */
    public final List<Residue> sequence;
    
    public final EnergyBreakdown energyBreakdown;

    /**
     * Public constructor.  To create a peptide, use the static factory method
     * createPeptide() instead.
     */
    public Peptide(String name, List<Atom> contents, SimpleWeightedGraph<Atom,DefaultWeightedEdge> connectivity,
                    List<Residue> sequence, EnergyBreakdown energyBreakdown)
    {
        super(name,contents,connectivity);
        this.sequence = sequence;
        this.energyBreakdown = energyBreakdown;
    }

    /**
     * Replaces one of the residues with a new one with the same atoms, but of a new type.
     * Meant for cis/trans proline changes.
     */
    public Peptide replaceResidueType(Residue oldResidue, Residue newResidue)
    {
        List<Residue> newSequence = new LinkedList<>();
        for (Residue r : sequence)
            {
                if ( r.equals(oldResidue) )
                    newSequence.add(newResidue);
                else
                    newSequence.add(r);
            }
        return new Peptide(name, contents, connectivity, newSequence, energyBreakdown);
    }
    
    /**
     * Factory method to create a peptide given a map of old atoms to new atoms.  Should be used
     * to move atoms.
     * @param atomMap a map from old atoms to new atoms (does not have to include all atoms)
     */
    public Peptide moveAtoms2(Map<Atom,Atom> atomMap)
    {
        // copy the list of vertices
        List<Atom> newContents = new LinkedList<Atom>();
        for (Atom a : contents)
            {
                if ( atomMap.containsKey(a) )
                    newContents.add(atomMap.get(a));
                else
                    newContents.add(a);
            }

        // populate a new connectivity graph
        SimpleWeightedGraph<Atom,DefaultWeightedEdge> newConnectivity = new SimpleWeightedGraph<Atom,DefaultWeightedEdge>(DefaultWeightedEdge.class);
        for (Atom newAtom : newContents)
            newConnectivity.addVertex(newAtom);
        for (DefaultWeightedEdge e : connectivity.edgeSet())
            {
                // get old edge data
                Double bondOrder = connectivity.getEdgeWeight(e);
                Atom fromAtom    = connectivity.getEdgeSource(e);
                Atom toAtom      = connectivity.getEdgeTarget(e);

                // replace any changes
                if ( atomMap.containsKey(fromAtom) )
                    fromAtom = atomMap.get(fromAtom);
                if ( atomMap.containsKey(toAtom)   )
                    toAtom   = atomMap.get(toAtom);

                // create new edge
                DefaultWeightedEdge newEdge = newConnectivity.addEdge(fromAtom,toAtom);
                newConnectivity.setEdgeWeight(newEdge, bondOrder);
            }

        // create new sequence/residues
        List<Residue> newSequence = new LinkedList<>();
        for (Residue r : sequence)
            newSequence.add(r.moveAtoms(atomMap));
        newSequence = ImmutableList.copyOf(newSequence);
        
        // return result
        // assumes the new peptide has the same name and won't have an EnergyBreakdown yet
        return new Peptide(name, newContents, newConnectivity, newSequence, EnergyBreakdown.BLANK);
    }

    /** returns a new peptide that is exactly the same but has a new name */
    public Peptide setName(String newName)
    {
        return new Peptide(newName, contents, connectivity, sequence, energyBreakdown);
    }

    /** returns a new peptide that has its EnergyBreakdown set */
    public Peptide setEnergyBreakdown(EnergyBreakdown energyBreakdown)
    {
        return new Peptide(name, contents, connectivity, sequence, energyBreakdown);
    }

    public Peptide setMolecule(Molecule newMolecule)
    {
        Molecule oldMolecule = (Molecule)this;
        Map<Atom,Atom> atomMap = getAtomMap(oldMolecule,newMolecule);
        return moveAtoms2(atomMap);
    }

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

    /** assumes that molecule1 and molecule2 are different conformations
     * returns a map from atoms in molecule 1 to atoms in molecule 2
     * map will not include atoms that don't move
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
     * Shakes the torsions until reasonable.
     * Won't shake hairpin residues.
     * @return the pseudo-minimized peptide
     */
    public Peptide shake()
    {
        ThreadLocalRandom random = ThreadLocalRandom.current();
        int numberOfResidues = sequence.size();

        // shake torsions until reasonable
        Peptide peptide = this;
        double energy = peptide.getOPLSenergy();

        for (int j=0; j<100; j++)
            {
                //System.out.println(j);

                // choose a random residue to mutate
                int residueIndex = random.nextInt(numberOfResidues);
                Residue residue = peptide.sequence.get(residueIndex);

                // don't shake any hairpin residues
                if ( residue.description.indexOf("hairpin") > -1 )
                    {
                        //System.out.println("skipping " + residue.description);
                        j--;
                        continue;
                    }

                // make mutations
                Peptide tempPeptide = BackboneMutator.mutatePhiPsi(peptide, residue);

                residue = tempPeptide.sequence.get(residueIndex);
                tempPeptide = BackboneMutator.mutateOmega(tempPeptide, residue);

                AminoAcid.RotamerType rotamerType = residue.aminoAcid.rotamerType;
                if ( rotamerType == AminoAcid.RotamerType.IS_ROTAMERIC ||
                      rotamerType == AminoAcid.RotamerType.NON_ROTAMERIC   )
                    {
                        residue = tempPeptide.sequence.get(residueIndex);
                        tempPeptide = RotamerMutator.mutateChis(tempPeptide, residue);
                    }

                if ( tempPeptide.checkCloseContacts() == false )
                    {
                        peptide = tempPeptide;
                        break;
                    }
                else
                    {
                        double thisEnergy = tempPeptide.getOPLSenergy();
                        if (thisEnergy < energy)
                            {
                                peptide = tempPeptide;
                                energy = thisEnergy;
                            }
                    }
            }

        return peptide;
    }

    /**
     * Compares peptides on the basis of their total energy in
     * energy breakdown (allows for sorting lists in ascending order) 
     * @param p2 the peptide to compare this peptide to
     */
    @Override
    public int compareTo(Peptide p2)
    {
        if ( this.energyBreakdown == null || p2.energyBreakdown == null )
            throw new NullPointerException("null energy breakdown -- cannot compare");
        return energyBreakdown.totalEnergy > p2.energyBreakdown.totalEnergy ? 1 : (energyBreakdown.totalEnergy < p2.energyBreakdown.totalEnergy ? -1 : 0);
    }

    public static void report(Peptide peptide, int j)
    {
        Residue r = peptide.sequence.get(j-1);
        System.out.println(r.toString(peptide));
    }

    @Override
    public int hashCode()
    {
        return Objects.hash(contents);
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

        Peptide anotherPeptide = (Peptide)obj;
        if (contents.equals(anotherPeptide.contents))
            return true;
       /* if ( contents.size() == anotherPeptide.contents.size() )
            {
                for ( int i=0; i < contents.size(); i++ )
                    {
                        Atom a1 = contents.get(i);
                        Atom a2 = anotherPeptide.contents.get(i);
                        if ( ! a1.equals(a2) )
                            System.out.println("mismatch at index " + i + ": " + a1.toString() + " / " + a2.toString());
                    }
                Scanner scanner = new Scanner(System.in);
                scanner.nextLine();
            }*/
        return false;
    }

    /** makes a bunch of peptides and checks their invariants */
    public static void main(String[] args)
    {
        System.out.println(OmegaLibrary.INSTANCE);
        System.out.println(RamachandranLibrary.INSTANCE);
        System.out.println(RotamerLibrary.getDescription());

        List<Integer> indexList = ImmutableList.of(1,2,5,6);
        List<Double> timings = new LinkedList<>();
        for ( int i=0; i < 100; i++ )
            {
                long startTime = System.currentTimeMillis();
                List<String> inputSequence = new LinkedList<>();
                for (int j=0; j<5; j++)
                    inputSequence.add(AminoAcid.getRandom());
                inputSequence.add("Dpro");
                inputSequence.add("Gly");
                for (int j=0; j<5; j++)
                    inputSequence.add(AminoAcid.getRandom());
                System.out.println(inputSequence);
                String[] array = inputSequence.toArray(new String[inputSequence.size()]);
                List<ProtoAminoAcid> preSequence = ProtoAminoAcidLibrary.getSequence(array,false);
                
                Peptide peptide = createPeptide(preSequence,false);
                //for (Residue r : peptide.sequence )
                //    System.out.println(r.description);
                
                double energy = peptide.getOPLSenergy();
                boolean tooClose = peptide.checkCloseContacts();

                for (int j=0; j < 100; j++)
                    {
                        System.out.print("iteration " + j + "\r");
                        Peptide lastPeptide = peptide;
                        
                        // select a residue at random to mutate
                        int residueIndex = indexList.get( ThreadLocalRandom.current().nextInt(indexList.size()) ) - 1;
                        Residue residue = peptide.sequence.get(residueIndex); 
                        
                        // make mutations
                        Peptide tempPeptide = BackboneMutator.mutatePhiPsi(peptide, residue);
                        
                        residue = tempPeptide.sequence.get(residueIndex);
                        tempPeptide = BackboneMutator.mutateOmega(tempPeptide, residue);

                        AminoAcid.RotamerType rotamerType = residue.aminoAcid.rotamerType;
                        if ( rotamerType == AminoAcid.RotamerType.IS_ROTAMERIC ||
                             rotamerType == AminoAcid.RotamerType.NON_ROTAMERIC   )
                            {
                                residue = tempPeptide.sequence.get(residueIndex);
                                tempPeptide = RotamerMutator.mutateChis(tempPeptide, residue);
                            }

                        if ( j > 20 && tempPeptide.checkCloseContacts() == false )
                            {
                                peptide = tempPeptide;
                                energy = tempPeptide.getOPLSenergy();
                                break;
                            }
                        else
                            {
                                double thisEnergy = tempPeptide.getOPLSenergy();
                                if ( thisEnergy < energy )
                                    {
                                        peptide = tempPeptide;
                                        energy = thisEnergy;
                                    }
                            }
                    }
                System.out.println();
                
                GaussianInputFile gjf = new GaussianInputFile(peptide);
                gjf.write(String.format("test_peptides/peptide_%04d.gjf", i));

                TinkerXYZInputFile tinkerFile = new TinkerXYZInputFile(peptide);
                tinkerFile.write(String.format("test_peptides/peptide_%04d.xyz", i));

                long endTime = System.currentTimeMillis();
                double elapsedTime = (endTime-startTime)/1000.0;
                timings.add(elapsedTime);

                System.out.println(String.format("Peptide %03d: elapsed = %.3f, close = %b, energy = %.2f\n", i, elapsedTime, tooClose, energy));
            }
        
        double averageTime = 0.0;
        for (Double d : timings)
            averageTime += d;
        averageTime = averageTime / timings.size();
        System.out.println(String.format("\nAverage time = %.3f s", averageTime));
    }
}
