import java.util.*;
import java.io.*;
import com.google.common.collect.*;
import org.apache.commons.math3.geometry.euclidean.threed.*;
import org.jgrapht.*;
import org.jgrapht.graph.*;
import org.jgrapht.alg.*;
import java.util.concurrent.*;


public class PeptideFactory 
{
    /** This class is not instantiable. */
    private PeptideFactory()
    {
        throw new IllegalArgumentException("not instantiable");
    }


    /** convenience method for peptide creation with hairpin formation set to true */
    public static Peptide createPeptide(List<ProtoAminoAcid> inputSequence)
    {
        return createPeptide(inputSequence,true);
    }

    /**
     * Builds a peptide given a sequence of ProtoAminoAcids.
     * The input sequence should have exactly one D-Pro-Gly unit in the N to C direction.
     * @param inputSequence the requested sequence specified in the N to C direction
     * @param hairpin set to true if you want auto-generation of a hairpin
     * @return the Peptide that embodies the geometry and metadata of the sequence
     */
    public static Peptide createPeptide(List<ProtoAminoAcid> inputSequence, boolean hairpin)
    {
        // check validity
        if ( inputSequence.size() < 2 )
            throw new IllegalArgumentException("A peptide must have at least two residues in it.");
        // check that this amino acid contains exactly one D-Pro

        // create temporary lists for creating a new Residue
        List<AminoAcid>          aminoAcids           = new LinkedList<>();
        List<ResidueType>        residueTypes         = new LinkedList<>();
        List<ProtoTorsion>       omegas               = new LinkedList<>();
        List<ProtoTorsion>       phis                 = new LinkedList<>();
        List<ProtoTorsion>       psis                 = new LinkedList<>();
        List<List<ProtoTorsion>> chis                 = new LinkedList<>();
        List<ProtoTorsion>       XHtorsions           = new LinkedList<>();
        List<List<ProtoTorsion>> frozenTorsions       = new LinkedList<>();
        List<List<Atom>>         atoms                = new LinkedList<>();
        List<List<Atom>>         frozenAtoms          = new LinkedList<>();
        List<Atom>               backboneHNs          = new LinkedList<>();
        List<Atom>               imidazoleHNs         = new LinkedList<>();
        List<Atom>               imidazoleNs          = new LinkedList<>();
        List<List<Atom>>         otherHNatoms         = new LinkedList<>();
        List<Double>             referenceEnergies    = new LinkedList<>();
        List<String>             descriptions         = new LinkedList<>();
        List<Pair<Atom,Atom>>    prochiralConnections = new LinkedList<>(); 

        // for Molecule part of new Peptide
        String                   newName              = "";
        List<Atom>               newContents          = new LinkedList<>();
        SimpleWeightedGraph<Atom,DefaultWeightedEdge> newConnectivity
                                                      = new SimpleWeightedGraph<>(DefaultWeightedEdge.class);
        Set<Pair<Atom,Atom>>     newBonds             = new HashSet<>();

        List<Atom>               allFrozenAtoms       = new LinkedList<>();

        // keep track of connections so new bonds can be formed
        List<Pair<Atom,Atom>>    NStickyConnections   = new LinkedList<>();
        List<Pair<Atom,Atom>>    CStickyConnections   = new LinkedList<>();
        List<Pair<Atom,Atom>>    newConnections       = new LinkedList<>();
        List<Pair<Integer,Integer>> newIndexConnections = new LinkedList<>();

        // build peptide in N to C direction
        for (int i=0; i < inputSequence.size(); i++)
            {
                // get data for this amino acid
                ProtoAminoAcid p = inputSequence.get(i);
                
                AminoAcid                tempAminoAcid           = p.r.aminoAcid;
                ResidueType              tempResidueType         = p.r.residueType;
                ProtoTorsion             tempOmega               = p.r.omega;
                ProtoTorsion             tempPhi                 = p.r.phi;
                ProtoTorsion             tempPsi                 = p.r.psi;
                List<ProtoTorsion>       tempChis                = p.r.chis;
                ProtoTorsion             tempXHtorsion           = p.r.XHtorsion;
                List<ProtoTorsion>       tempFrozenTorsions      = p.r.frozenTorsions;
                List<Atom>               tempAtoms               = p.r.atoms;
                List<Atom>               tempFrozenAtoms         = p.r.frozenAtoms;
                Atom                     tempBackboneHN          = p.r.backboneHN;
                Atom                     tempImidazoleHN         = p.r.imidazoleHN;
                Atom                     tempImidazoleN          = p.r.imidazoleN;
                List<Atom>               tempOtherHNatoms        = p.r.otherHNatoms;
                Double                   tempReferenceEnergy     = p.r.referenceEnergy;
                String                   tempDescription         = p.r.description;
                Pair<Atom,Atom>          tempNStickyConnection   = p.NStickyConnection;
                Pair<Atom,Atom>          tempCStickyConnection   = p.CStickyConnection;
                Pair<Atom,Atom>          tempProchiralConnection = p.r.prochiralConnection;
                
                Molecule                 tempMolecule            = p.molecule;

                // append to name
                if ( i == 0 )
                    newName = newName + tempAminoAcid.shortName;
                else
                    newName = newName + " - " + tempAminoAcid.shortName;

                // add header data
                aminoAcids.add(tempAminoAcid);
                residueTypes.add(tempResidueType);

                // add relevant atoms from this ProtoAminoAcid
                tempAtoms = new LinkedList<>(tempAtoms);
                if ( i > 0 )
                    {
                        // delete N-terminal cap
                        Atom keepAtom = tempNStickyConnection.getFirst();
                        Atom deleteAtom = tempNStickyConnection.getSecond();
                        Set<Atom> toBeDeleted = tempMolecule.getHalfGraph(keepAtom, deleteAtom);
                        tempAtoms.removeAll(toBeDeleted);
                    }

                if ( i < inputSequence.size() - 1 )
                    {
                        // delete C-terminal cap
                        Atom keepAtom = tempCStickyConnection.getFirst();
                        Atom deleteAtom = tempCStickyConnection.getSecond();
                        Set<Atom> toBeDeleted = tempMolecule.getHalfGraph(keepAtom, deleteAtom);
                        tempAtoms.removeAll(toBeDeleted);
                    }
                newContents.addAll(tempAtoms);
                atoms.add(tempAtoms);

                // add relevant connectivity from each ProtoAminoAcid
                for (DefaultWeightedEdge e : tempMolecule.connectivity.edgeSet())
                    {
                        // assumes all edge weights are 1
                        Atom fromAtom = tempMolecule.connectivity.getEdgeSource(e);
                        Atom toAtom   = tempMolecule.connectivity.getEdgeTarget(e);
                        if ( tempAtoms.contains(fromAtom) && tempAtoms.contains(toAtom) )
                            newBonds.add( new Pair<Atom,Atom>(fromAtom, toAtom) );
                    }

                // form new bond
                NStickyConnections.add(tempNStickyConnection);
                CStickyConnections.add(tempCStickyConnection);

                if ( i > 0 )  // no need to form a bond on the first residue
                    {
                        // join the C sticky atom of the previous residue
                        // to the N sticky atom of the current residue
                        // sticky atom connections are defined as ordered pairs: (keep, delete)
                        Atom lastCStickyAtom = CStickyConnections.get(i-1).getFirst();
                        Atom thisNStickyAtom = NStickyConnections.get(i).getFirst();
                        Pair<Atom,Atom> newConnection = new Pair<>(lastCStickyAtom, thisNStickyAtom);
                        //System.out.println(newContents.indexOf(lastCStickyAtom)+1);
                        //System.out.println(newContents.indexOf(thisNStickyAtom)+1);

                        // mark this in a list of all the new peptide bonds
                        newConnections.add(newConnection);

                        // mark for addition to connectivity graph
                        newBonds.add(newConnection);

                        // specify the new connections by atom number because atoms will be moved
                        // indices: 0, 1, 2, ..., N
                        Integer fromAtomIndex = newContents.indexOf(lastCStickyAtom);
                        Integer toAtomIndex = newContents.indexOf(thisNStickyAtom);
                        Pair<Integer,Integer> newIndexPair = new Pair<>(fromAtomIndex, toAtomIndex);
                        newIndexConnections.add(newIndexPair);
                    }

                // update backbone angles
                if ( i == 0 )
                    {
                        // copy ProtoTorsions for first residue
                        omegas.add(tempOmega);
                        phis.add(tempPhi);
                        psis.add(tempPsi);

                        // copy frozen torsions
                        frozenTorsions.add(tempFrozenTorsions);
                    }
                else if ( i > 0 )
                    {
                        // get last ProtoTorsions
                        ProtoTorsion lastPsi = psis.remove(psis.size()-1);
                        
                        // update ProtoTorsions for previous residue
                        ProtoTorsion newLastPsi = new ProtoTorsion(lastPsi.atom1, lastPsi.atom2, lastPsi.atom3, tempOmega.atom3);
                        psis.add(newLastPsi);

                        // update ProtoTorsions for current residue
                        ProtoTorsion newOmega = new ProtoTorsion(lastPsi.atom2, lastPsi.atom3, tempOmega.atom3, tempOmega.atom4);
                        ProtoTorsion newPhi   = new ProtoTorsion(lastPsi.atom3, tempPhi.atom2,   tempPhi.atom3,   tempPhi.atom4);
                        omegas.add(newOmega);
                        phis.add(newPhi);
                        psis.add(tempPsi);

                        // correct the last residue's frozen psi torsion if necessary
                        List<ProtoTorsion> lastFrozenTorsions = frozenTorsions.remove(frozenTorsions.size()-1);
                        List<ProtoTorsion> newLastFrozenTorsions = new LinkedList<>();
                        for (ProtoTorsion t : lastFrozenTorsions)
                            {
                                if ( t.equals(lastPsi) )
                                    newLastFrozenTorsions.add(newLastPsi);
                                else
                                    newLastFrozenTorsions.add(t);
                            }
                        frozenTorsions.add(newLastFrozenTorsions);

                        // adjust current residue's frozen torsions if necessary
                        List<ProtoTorsion> tempFrozenTorsions2 = new LinkedList<>();
                        for (ProtoTorsion t : tempFrozenTorsions)
                            {
                                if ( t.equals(p.r.phi) )
                                    {
                                        //System.out.println(tempAminoAcid.toString() + " phi " + t.toString(tempMolecule));;
                                        tempFrozenTorsions2.add(newPhi);
                                    }
                                else if ( t.equals(p.r.psi) )
                                    {
                                        //System.out.println(tempAminoAcid.toString() + " psi " + t.toString(tempMolecule));;
                                        tempFrozenTorsions2.add(tempPsi);
                                    }
                                else if ( t.equals(p.r.omega) )
                                    {
                                        //System.out.println(tempAminoAcid.toString() + " omega " + t.toString(tempMolecule));;
                                        tempFrozenTorsions2.add(newOmega);
                                    }
                                else
                                    {
                                        //System.out.println(tempAminoAcid.toString() + " chi " + t.toString(tempMolecule));;
                                        tempFrozenTorsions2.add(t);
                                    }
                            }
                        
                        // update list of lists of frozen torsions
                        frozenTorsions.add(tempFrozenTorsions2);
                    }

                // add sidechain torsion angles
                chis.add(tempChis);
                XHtorsions.add(tempXHtorsion);

                // add other fields
                frozenAtoms.add(tempFrozenAtoms);
                backboneHNs.add(tempBackboneHN);
                imidazoleHNs.add(tempImidazoleHN);
                imidazoleNs.add(tempImidazoleN);
                otherHNatoms.add(tempOtherHNatoms);
                referenceEnergies.add(tempReferenceEnergy);
                descriptions.add(tempDescription);
                prochiralConnections.add(tempProchiralConnection);
            }

        // create new connectivity graph
        for (Atom a : newContents)
            newConnectivity.addVertex(a);
        for (Pair<Atom,Atom> p : newBonds)
            {
                // assumes all edges have a weight of 1.0
                Atom fromAtom = p.getFirst();
                Atom toAtom = p.getSecond();
                DefaultWeightedEdge e = newConnectivity.addEdge(fromAtom, toAtom);
                newConnectivity.setEdgeWeight(e, 1.0);
                //System.out.println( (newContents.indexOf(fromAtom)+1) + " -- " + (newContents.indexOf(toAtom)+1) );
            }

        // create residue objects to make a Sequence
        List<Residue> sequence = new LinkedList<>();
        for (int i=0; i < inputSequence.size(); i++)
            {
                Residue residue = new Residue(aminoAcids.get(i),
                                              residueTypes.get(i),
                                              omegas.get(i),
                                              phis.get(i),
                                              psis.get(i),
                                              chis.get(i),
                                              XHtorsions.get(i),
                                              frozenTorsions.get(i),
                                              atoms.get(i),
                                              frozenAtoms.get(i),
                                              backboneHNs.get(i),
                                              imidazoleHNs.get(i),
                                              imidazoleNs.get(i),
                                              otherHNatoms.get(i),
                                              referenceEnergies.get(i),
                                              descriptions.get(i),
                                              prochiralConnections.get(i),
                                              inputSequence.get(i));

                sequence.add(residue);
            }
        sequence = ImmutableList.copyOf(sequence);

        // create new peptide object
        // no constraints / energyBreakdown by default
        EnergyBreakdown energyBreakdown = EnergyBreakdown.BLANK;
        List<DistanceConstraint> constraints = ImmutableList.of();
        Peptide newPeptide = new Peptide(newName, newContents, newConnectivity, sequence, energyBreakdown, allFrozenAtoms, constraints);

        // adjust distances
        Peptide movedPeptide = newPeptide;
        for (Pair<Integer,Integer> bond : newIndexConnections)
            {
                Atom fromAtom = movedPeptide.contents.get(bond.getFirst());
                Atom toAtom   = movedPeptide.contents.get(bond.getSecond());
                movedPeptide  = movedPeptide.setMolecule( movedPeptide.setDistance(fromAtom,toAtom,1.32) );
            }
        
        // set peptide bonds to sp2
        for (int i=1; i < inputSequence.size(); i++)
            {
                // set amide carbonyl carbon to sp2
                Residue r = movedPeptide.sequence.get(i);
                ProtoTorsion omega = r.omega;
                Atom fromAtom = omega.atom2;
                Atom toAtom = omega.atom3;
                movedPeptide = movedPeptide.setMolecule( movedPeptide.set_sp2(fromAtom,toAtom) );

                // set amide nitrogen to sp2
                r = movedPeptide.sequence.get(i);
                omega = r.omega;
                fromAtom = omega.atom3;
                toAtom = omega.atom2;

                // only force (amide substituent 1 -- N -- amide substituent 2) angle to 120
                // if this is not a proline
                if ( r.aminoAcid.isProline() )
                    movedPeptide = movedPeptide.setMolecule( movedPeptide.set_sp2(fromAtom,toAtom,false) );
                else
                    movedPeptide = movedPeptide.setMolecule( movedPeptide.set_sp2(fromAtom,toAtom) );
            }

        // set backbone torsions to trans and set hairpin
        for (int i=0; i < inputSequence.size(); i++)
            {
                Residue r = movedPeptide.sequence.get(i);
                ProtoTorsion omega = r.omega;
                
                if ( hairpin && r.aminoAcid == AminoAcid.DPRO )
                    {
                        // fix hairpin geometry
                        movedPeptide = movedPeptide.setMolecule( movedPeptide.setDihedral(omega, -176.0) );

                        r = movedPeptide.sequence.get(i);
                        ProtoTorsion psi = r.psi;
                        movedPeptide = movedPeptide.setMolecule( movedPeptide.setDihedral(psi,-132.0) );

                        i++;
                        r = movedPeptide.sequence.get(i);
                        if ( r.aminoAcid != AminoAcid.GLY )
                            throw new IllegalArgumentException("glycine is expected to follow D-proline");
                        omega = r.omega;
                        movedPeptide = movedPeptide.setMolecule( movedPeptide.setDihedral(omega,175.0) );

                        r = movedPeptide.sequence.get(i);
                        ProtoTorsion phi = r.phi;
                        movedPeptide = movedPeptide.setMolecule( movedPeptide.setDihedral(phi,-96.0) );

                        r = movedPeptide.sequence.get(i);
                        psi = r.psi;
                        movedPeptide = movedPeptide.setMolecule( movedPeptide.setDihedral(psi,9.0) );
                    }
                else
                    // set to trans peptide bonds
                    movedPeptide = movedPeptide.setMolecule( movedPeptide.setDihedral(omega, 180.0) );
            }

        // set proline type for any prolines
        for (int i=0; i < inputSequence.size(); i++)
            {
                Residue oldResidue = movedPeptide.sequence.get(i);
                if ( oldResidue.aminoAcid == AminoAcid.PRO )
                    {
                        Residue newResidue = oldResidue.setProline();
                        movedPeptide = movedPeptide.replaceResidueType(oldResidue,newResidue);
                    }
            }

        // return final result
        Peptide finalPeptide = movedPeptide;
        return finalPeptide;
    }

    /**
     * Computes a rough Lennard-Jones steric energy for this molecule.  Answer is normalized
     * by the number of atoms; i.e. kcal/mol divided by the number of atoms.  Atoms that are
     * separated by one or two bonds are ignored.  Atoms separated by more than Settings.CUTOFF_DISTANCE
     * are also ignored.
     * @return the steric energy
     */
    public static double getOPLSenergy(Peptide peptide)
    {
        double energy = 0.0;
        for (int i=0; i < contents.size(); i++)
            {
                Atom atom1 = contents.get(i);
                Vector3D atom1position = atom1.position;
                for (int j=i+1; j < contents.size(); j++)
                    {
                        Atom atom2 = contents.get(j);
                        // ignore if atoms are too close in the connectivity graph
                        if ( ! areSeparated(atom1,atom2) )
                            continue;

                        // ignore if atoms are too far apart
                        Vector3D atom2position = atom2.position;
                        double distance = Vector3D.distance(atom1position, atom2position);
                        if ( distance > Settings.CUTOFF_DISTANCE )
                            continue;

                        // prevent overflow
                        if ( distance < 0.5 )
                            distance = 0.5;

                        // get parameters
                        double epsilon1 = atom1.element.epsilon;
                        double sigma1   = atom1.element.sigma;
                        double epsilon2 = atom2.element.epsilon;
                        double sigma2   = atom2.element.sigma;

                        // apply combination rules
                        double epsilon = epsilon1;
                        double sigma   = sigma1;
                        if ( epsilon1 != epsilon2 )
                            epsilon = Math.sqrt(epsilon1 * epsilon2);
                        if ( sigma1 != sigma2 )
                            sigma = Math.sqrt(sigma1 * sigma2);

                        // compute energy
                        double temp = Math.pow(sigma/distance, 6);
                        energy += 4.0 * epsilon * temp * (temp - 1.0);
                    }
            }
        return energy / contents.size();
    }


    /**
     * Creates a peptide with a non-frozen hairpin for reference energy calculations.
     * Torsions are shaken so that there are no clashes.
     * @param numberOfResidues the number of residues on each arm, such that the total
     *        number of residues will be 2*numberOfResidues+2
     * @return the random well-formed peptide
     */
    public static Peptide createReferencePeptide(int numberOfResidues)
    {
        // calculate peptide length
        if ( numberOfResidues < 1 )
            throw new IllegalArgumentException("number of residues in each arm must be greater than 1");
        int peptideLength = 2*numberOfResidues + 2;

        // generate a random sequence
        List<String> inputSequence = new LinkedList<>();
        for (int i=0; i < numberOfResidues; i++)
            inputSequence.add(AminoAcid.getRandom());
        inputSequence.add("dpro");
        inputSequence.add("gly");
        for (int i=0; i < numberOfResidues; i++)
            inputSequence.add(AminoAcid.getRandom());
        
        String[] array = inputSequence.toArray(new String[inputSequence.size()]);
        List<ProtoAminoAcid> preSequence = ProtoAminoAcidLibrary.getSequence(array,false);
        System.out.println();
        for (ProtoAminoAcid p : preSequence)
            System.out.print(p.r.description + ", ");
        System.out.println("\n");

        Peptide peptide = Peptide.createPeptide(preSequence,false);
        double energy = peptide.getOPLSenergy();

        // shake torsions until reasonable
        ThreadLocalRandom random = ThreadLocalRandom.current();
        for (int j=0; j<100; j++)
            {
                Peptide lastPeptide = peptide;
                
                // choose a random residue to mutate
                int residueIndex = random.nextInt(numberOfResidues);
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

                if ( tempPeptide.checkCloseContacts() == false )
                    {
                        peptide = tempPeptide;
                        energy = tempPeptide.getOPLSenergy();
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
                
                Peptide peptide = PeptideFactory.createPeptide(preSequence,false);
              
                GaussianInputFile gjf = new GaussianInputFile(peptide);
                gjf.write(String.format("test_peptides/peptide_%04d.gjf", i));

                long endTime = System.currentTimeMillis();
                double elapsedTime = (endTime-startTime)/1000.0;
                timings.add(elapsedTime);

                System.out.println(String.format("Peptide %03d: elapsed = %.3f\n", i, elapsedTime));
            }
    }
}  
