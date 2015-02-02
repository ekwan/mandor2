import org.apache.commons.math3.geometry.euclidean.threed.*;
import java.io.*;
import java.util.*;
import java.util.concurrent.*;
import com.google.common.collect.*;
import org.jgrapht.*;
import org.jgrapht.graph.*;
import org.jgrapht.alg.*;

/**
 * This collects together methods for taking a peptide and mutating one of
 * its residues from one amino acid to another.
 *
 * Note that this class cannot currently convert between D and L amino acids.
 * The chirality of the amino acid is preserved at present because the old
 * sidechain is simply replaced with the new sidechain.  I don't expect this
 * kind of change, so I'm just going to throw an exception if this happens.
 * A future implementation will have to deal with glycine and other achiral
 * versions as well.
 */
public class SidechainMutator implements Mutator
{
    /** This class is not instantiable. */
    private SidechainMutator()
    {
        throw new IllegalArgumentException("not instantiable");
    }

    /**
     * Mutates the a residue in a given peptide to a specified amino acid.  No backbone atoms will be moved as a result of
     * a call to mutateSidechain, except in the case of moving from a non-proline to a proline or vice versa.  When that
     * happens, the HN will be altered, since non-prolines have one and prolines don't.  However, none of the other backbone
     * atoms will ever be moved.  After the mutation process, a rotamer is randomly chosen given the backbone context.
     * @param inputPeptide the peptide to mutate
     * @param inputResidue a residue inside inputPeptide that should be changed
     * @param targetPAA the template for the amino acid we want to mutate inputResidue into
     * @return the mutated peptide
     */
    public static Peptide mutateSidechain(Peptide inputPeptide, Residue inputResidue, ProtoAminoAcid targetPAA)
    {
        // check that inputResidue belongs to inputPeptide
        if ( ! inputPeptide.sequence.contains(inputResidue) )
            throw new IllegalArgumentException("inputResidue is not contained in inputPeptide");

        // if the residue is already the correct template, then just return the input peptide
        if ( inputResidue.description.equals(targetPAA.residue.description) )
            return inputPeptide;

        // don't allow changing between chiralities
        if ( ( inputResidue.aminoAcid.chirality == Chirality.L && targetPAA.residue.aminoAcid.chirality == Chirality.D ) ||
             ( inputResidue.aminoAcid.chirality == Chirality.D && targetPAA.residue.aminoAcid.chirality == Chirality.L )    )
            throw new IllegalArgumentException("chirality interchange forbidden -- see javadoc for class");
        if ( inputResidue.aminoAcid.chirality == Chirality.ACHIRAL && targetPAA.residue.aminoAcid.chirality == Chirality.D )
            throw new IllegalArgumentException("can't handle mutations to D yet");

        // this is the index of the residue we will be mutating
        int residueIndex = inputPeptide.sequence.indexOf(inputResidue);

        // set the backbone atoms to the types in targetPAA, but don't change their positions
        //
        // note that at this point, we have created a new peptide, so it is safe to modify the connectivity graph within
        // the scope of this method
        Peptide newPeptide = BackboneAnalysis.adjustBackboneAtomTypes(inputPeptide, inputResidue, targetPAA);

        // if switching from proline to non-proline, disconnect the sidechain-backbone connection and add the HN
        if ( inputResidue.aminoAcid.isProline() && !targetPAA.residue.aminoAcid.isProline() )
            newPeptide = convertProlineToNonProline(newPeptide, residueIndex, targetPAA);

        // if switching from non-proline to proline, delete the current HN
        else if ( !inputResidue.aminoAcid.isProline() && targetPAA.residue.aminoAcid.isProline() )
            newPeptide = convertNonProlineToProline(newPeptide, residueIndex);

        // special case for switching from proline to proline, only make disconnection between sidechain and backbone  
        else if ( inputResidue.aminoAcid.isProline() && targetPAA.residue.aminoAcid.isProline() )
            newPeptide = convertProlineToProline(newPeptide, residueIndex);

        // get the new sidechain atoms
        Sidechain newSidechain = Sidechain.makeNewSidechain(newPeptide, residueIndex, targetPAA);
        
        // delete the current sidechain atoms
        //System.out.println("Source residue: " + inputResidue);
        //System.out.println("Target PAA: " + targetPAA.r.aminoAcid);
        newPeptide = deleteCurrentSidechain(newPeptide, residueIndex);

        // add the new sidechain
        newPeptide = addNewSidechain(newPeptide, residueIndex, targetPAA, newSidechain);

        // create new peptide name
        String newName = generateName(newPeptide);

        // choose an appropriate rotamer for the new residue
        newPeptide = RotamerMutator.mutateChis(newPeptide, residueIndex);

        // create and return new peptide
        return new Peptide(newName, newPeptide.contents, newPeptide.connectivity, newPeptide.sequence, EnergyBreakdown.BLANK);
    }

    /**
     * Converts a proline to a non-proline.  Deletes the N-CD connection and adds the HN.
     * @return the modified peptide
     */
    private static Peptide convertProlineToNonProline(Peptide peptide, int residueIndex, ProtoAminoAcid targetPAA)
    {
        // get fields
        Residue residue = peptide.sequence.get(residueIndex);   // the residue we need to mutate
        Atom atomCD = residue.chis.get(1).atom4;                // the delta carbon
        Vector3D atomCDposition = atomCD.position;              // location of the delta carbon
        Atom atomN = residue.chis.get(0).atom1;                 // the peptide nitrogen 

        // figure out where to put the new HN
        double desiredDistance = 1.02;                          // N-H bond is 1.02 A
        double currentDistance = Vector3D.distance(atomN.position, atomCDposition);
        Vector3D translateVector = atomCDposition.subtract(atomN.position);
        double scaling = (desiredDistance - currentDistance) / currentDistance;

        Vector3D requiredTranslation = translateVector.scalarMultiply(scaling);
        Vector3D newHNPosition = atomCDposition.add(requiredTranslation);
        
        // create the HN atom
        Atom templateHN = targetPAA.residue.HN;
        Atom atomHN = new Atom(templateHN.element, newHNPosition, templateHN.type1, templateHN.type2, templateHN.surfaceTension);

        // add HN to the list of atoms
        List<Atom> newContents = new ArrayList<>(peptide.contents);
        newContents.add(atomHN);
        
        // add HN to a new residue
        List<Atom> updatedAtoms = new ArrayList<>(residue.atoms);
        updatedAtoms.add(atomHN);
        Residue newResidue = new Residue(residue.aminoAcid, residue.omega, residue.phi, residue.psi, residue.chis,
                                         atomHN, residue.N, residue.O, residue.C, residue.CA, residue.HA, residue.description,
                                         residue.prochiralConnection, updatedAtoms, residue.isHairpin);

        // replace the residue in the sequence
        List<Residue> newResidueList = new ArrayList<>(peptide.sequence);
        newResidueList.set(residueIndex,newResidue);
        
        // adjust the connectivity to add the N-HN bond --> note that this modifies the backing peptide
        // this is safe, because the backing peptide is meant to be a hidden intermediate local to this class
        SimpleWeightedGraph<Atom,DefaultWeightedEdge> newConnectivity = peptide.connectivity;
        newConnectivity.addVertex(atomHN);
        newConnectivity.addEdge(atomN, atomHN);
        
        // remove the N-CD connection
        DefaultWeightedEdge oldEdge = newConnectivity.removeEdge(atomN, atomCD);
        if ( oldEdge == null )
            throw new NullPointerException("failed to remove N-CD connection");

        // return the new peptide
        return new Peptide(peptide.name, newContents, newConnectivity, newResidueList, EnergyBreakdown.BLANK);
    }

    /**
     * Converts a non-proline to a proline.  Removes the HN.  The returned peptide will use the same
     * connectivity graph as the original.  That is, the original connectivity will be modified.
     */
    private static Peptide convertNonProlineToProline(Peptide peptide, int residueIndex)
    {
        // get fields
        Residue residue = peptide.sequence.get(residueIndex);
        Atom atomHN = residue.HN;
        if ( atomHN == null )
            throw new IllegalArgumentException("HN not found");

        // create new list of atoms
        List<Atom> contents = new ArrayList<>(peptide.contents);
        boolean success = contents.remove(atomHN);
        if ( !success )
            throw new IllegalArgumentException("failed to remove HN");

        // create new residue
        List<Atom> updatedAtoms = new ArrayList<>(residue.atoms);
        success = updatedAtoms.remove(atomHN);
        if ( !success)
            throw new IllegalArgumentException("failed to remove HN");
        Residue newResidue = new Residue(residue.aminoAcid, residue.omega, residue.phi, residue.psi, residue.chis,
                                         null, residue.N, residue.O, residue.C, residue.CA, residue.HA, residue.description,
                                         residue.prochiralConnection, updatedAtoms, residue.isHairpin);
 
        // adjust sequence
        List<Residue> newSequence = new ArrayList<>(peptide.sequence);
        newSequence.set(residueIndex,newResidue);

        // adjust connectivity --> note that this modifies the peptide
        SimpleWeightedGraph<Atom,DefaultWeightedEdge> connectivity = peptide.connectivity;
        connectivity.removeVertex(atomHN);

        // make the new peptide
        return new Peptide(peptide.name, contents, connectivity, newSequence, EnergyBreakdown.BLANK);
    }

    /**
    * Converts a proline to a proline. The returned peptide will have its ring opened to allow for removal of sidechain.
    * @param peptide the peptide that will have its proline residue modified
    * @param residueIndex the residue index of a proline in the peptide
    * @return the peptide with the proline ring opened
    */
    private static Peptide convertProlineToProline(Peptide peptide, int residueIndex)
    {
        //remove the C-ND connection
        Residue residue = peptide.sequence.get(residueIndex);
        Atom atomCD = residue.chis.get(1).atom4;
        Atom atomN = residue.chis.get(0).atom1;

        // change connectivity table -- note modifies peptide
        SimpleWeightedGraph<Atom,DefaultWeightedEdge> newConnectivity = peptide.connectivity;        
        DefaultWeightedEdge oldEdge = newConnectivity.removeEdge(atomN, atomCD);
        if ( oldEdge == null )
            throw new NullPointerException("failed to remove N-CD connection");
        
        //residue stays the same
        return new Peptide(peptide.name, peptide.contents, newConnectivity, peptide.sequence, EnergyBreakdown.BLANK);
    }

    /**
     * Creates a new name for a peptide.
     * @param peptide the peptide to analyze
     * @return the new name
     */
    public static String generateName(Peptide peptide)
    {
        String name = String.format("%5s", peptide.sequence.get(0).aminoAcid.shortName);
        for (int i=1; i < peptide.sequence.size(); i++)
            name += String.format(" - %5s", peptide.sequence.get(i).aminoAcid.shortName);
        return name;
    }

    /**
     * Deletes the sidechain of the specified peptide.  This will leave the specified residue in an inconsistent state
     * because it won't have any of its sidechain atoms.  This will also modify the peptide connectivity graph object
     * directly.
     * @return the modified peptide
     */
    private static Peptide deleteCurrentSidechain(Peptide peptide, int residueIndex)
    {
        // get the current sidechain atoms
        // assumes there's no issue with proline
        Residue residue = peptide.sequence.get(residueIndex);
        Pair<Atom,Atom> prochiralConnection = residue.prochiralConnection;
        //System.out.println(peptide.sequence.get(residueIndex));
        Set<Atom> sidechainAtoms = peptide.getHalfGraph(prochiralConnection.getFirst(), prochiralConnection.getSecond());

        // delete the sidechain atoms from the atom list
        List<Atom> contents = new ArrayList<>(peptide.contents);
        for (Atom a : sidechainAtoms)
            contents.remove(a);

        // update the connectivity table
        SimpleWeightedGraph<Atom, DefaultWeightedEdge> connectivity = peptide.connectivity;
        for (Atom a : sidechainAtoms)
            {
                boolean success = connectivity.removeVertex(a);
                if ( !success )
                    throw new IllegalArgumentException("failed to remove atom from connectivity when deleting sidechain");
            }

        // update the residue
        List<Atom> updatedAtoms = new ArrayList<>(residue.atoms);
        for (Atom a : sidechainAtoms)
            {
                boolean success = updatedAtoms.remove(a);
                if ( !success )
                    throw new IllegalArgumentException("failed to remove atom from residue when deleting sidechain");
            }
        List<ProtoTorsion> newChis = new ArrayList<>();
        Pair<Atom,Atom> dummyProchiralConnection = new Pair<>(residue.CA, residue.CA);
        Residue newResidue = new Residue(residue.aminoAcid, residue.omega, residue.phi, residue.psi, newChis,
                                         residue.HN, residue.N, residue.O, residue.C, residue.CA, residue.HA, residue.description,
                                         dummyProchiralConnection, updatedAtoms, residue.isHairpin);
 
        // update the sequence
        List<Residue> newSequence = new ArrayList<>(peptide.sequence);
        newSequence.set(residueIndex, newResidue);

        // make the new peptide
        return new Peptide(peptide.name, contents, connectivity, newSequence, EnergyBreakdown.BLANK);
    }

    /**
     * Adds the specified sidechain to the given peptide.  Adds the sidechain atoms to the contents and connectivity
     * graph.  Also adds the needed connections between the backbone and the sidechain.  Assumes the old sidechain has
     * already been deleted from the peptide.  Will modify peptide.connectivity directly.
     * @param peptide the peptide to add the sidechain to
     * @param residueIndex the sequence position at which to add the sidechain
     * @param targetPAA what we want to turn the residue into
     * @param sidechain the new sidechain
     * @return the modified peptide
     */
    private static Peptide addNewSidechain(Peptide peptide, int residueIndex, ProtoAminoAcid targetPAA, Sidechain sidechain)
    {
        // add the sidechain atoms to the atom list
        Residue residue = peptide.sequence.get(residueIndex);
        List<Atom> contents = new ArrayList<>(peptide.contents);
        for (Atom a : sidechain.contents)
            contents.add(a);

        // add the sidechain atoms to the connectivity
        SimpleWeightedGraph<Atom,DefaultWeightedEdge> connectivity = peptide.connectivity;
        for (Atom a : sidechain.contents)
            {
                boolean success = connectivity.addVertex(a);
                if ( !success )
                    throw new IllegalArgumentException("failed to add a new sidechain atom to the peptide connectivity because it is already there");
            }

        // form the desired connections between the backbone and the sidechain
        connectivity.addEdge(residue.CA, sidechain.atomCB);
        if ( sidechain.atomCD != null )
            connectivity.addEdge(residue.N, sidechain.atomCD);

        // add the connections internal to the sidechain
        for (DefaultWeightedEdge e : sidechain.connectivity.edgeSet())
            {
                Atom fromAtom = connectivity.getEdgeSource(e);
                Atom toAtom   = connectivity.getEdgeTarget(e);
                connectivity.addEdge(fromAtom,toAtom);
            }

        // make new chis
        List<ProtoTorsion> chis = new ArrayList<>();               // the new chis
        List<Atom> chiAtoms = sidechain.chiAtoms;                  // CB, CG, CD, ... in the new sidechain
        
        if ( targetPAA.residue.aminoAcid.isProline() )
            {
                ProtoTorsion newChi1 = new ProtoTorsion(residue.N, residue.CA, chiAtoms.get(0), chiAtoms.get(1));
                ProtoTorsion newChi2 = new ProtoTorsion(residue.CA, chiAtoms.get(0), chiAtoms.get(1), chiAtoms.get(2));
                ProtoTorsion newChi3 = new ProtoTorsion(chiAtoms.get(0), chiAtoms.get(1), chiAtoms.get(2), residue.N);
                chis.add(newChi1);
                chis.add(newChi2);
                chis.add(newChi3);
            }
        else
            {
                if ( targetPAA.residue.chis.size() >= 1 ) // deal with chi1
                    {
                        // we need the old N and CA, but the new CB and CG
                        ProtoTorsion newChi1 = new ProtoTorsion(residue.N, residue.CA, chiAtoms.get(0), chiAtoms.get(1));
                        chis.add(newChi1);
                    }

                if ( targetPAA.residue.chis.size() >= 2 ) // deal with chi2
                    {
                        // we need the old N and CA, but the new CB and CG
                        ProtoTorsion newChi2 = new ProtoTorsion(residue.CA, chiAtoms.get(0), chiAtoms.get(1), chiAtoms.get(2));
                        chis.add(newChi2);
                    }

                if ( targetPAA.residue.chis.size() >= 3 ) // deal with all higher order chis
                    {
                        // these are only made of the new atoms
                        for (int i=2; i < targetPAA.residue.chis.size(); i++)
                            {
                                ProtoTorsion newChiN = new ProtoTorsion(chiAtoms.get(i-2), chiAtoms.get(i-1), chiAtoms.get(i), chiAtoms.get(i+1));
                                chis.add(newChiN);
                            }
                    }
            }
        if ( chis.size() != targetPAA.residue.chis.size() )
            throw new IllegalArgumentException("chi size mismatch for " + targetPAA.residue.description);

        // get fields for a new residue
        List<Atom> residueAtoms = new ArrayList<>(residue.atoms);
        for (Atom a : sidechain.contents)                          // add the sidechain atoms to the residue atoms
            residueAtoms.add(a);
        Pair<Atom,Atom> newProchiralConnection = new Pair<>(residue.CA, sidechain.atomCB);

        // optional checks
        checkProtoTorsion(contents, residue.omega);
        checkProtoTorsion(contents, residue.phi);
        checkProtoTorsion(contents, residue.psi);
        for (ProtoTorsion torsion : chis)
            checkProtoTorsion(contents,torsion);
        for (Atom a : residueAtoms)
            if (! contents.contains(a) )
                throw new IllegalArgumentException("residue-contents mismatch");

        // make a new residue
        Residue newResidue = new Residue(targetPAA.residue.aminoAcid, residue.omega, residue.phi, residue.psi, chis,
                                         residue.HN, residue.N, residue.O, residue.C, residue.CA, residue.HA, targetPAA.residue.description,
                                         newProchiralConnection, residueAtoms, false);

        // update the sequence
        List<Residue> newSequence = new ArrayList<>(peptide.sequence);
        newSequence.set(residueIndex,newResidue);

        // make the new peptide
        return new Peptide(peptide.name, contents, connectivity, newSequence, EnergyBreakdown.BLANK);
    }

    /**
     * Ensures that the atoms in a ProtoTorsion in a some list of atoms.  Throws an exception if they aren't.
     */
    public static void checkProtoTorsion(List<Atom> contents, ProtoTorsion torsion)
    {
        if ( ! contents.contains(torsion.atom1) )
            throw new IllegalArgumentException("atom1 not in contents");
        else if ( ! contents.contains(torsion.atom2) )
            throw new IllegalArgumentException("atom2 not in contents");
        else if ( ! contents.contains(torsion.atom3) )
            throw new IllegalArgumentException("atom3 not in contents");
        else if ( ! contents.contains(torsion.atom4) )
            throw new IllegalArgumentException("atom4 not in contents");
    }

    /**
     * Represents the backbone atom assignments for a given residue.
     */
    public static class BackboneAnalysis implements Immutable
    {
        public final Atom atomHN, atomN, atomCA, atomC, atomO, atomHA;

        /**
         * This private constructor just copies fields.
         */
        private BackboneAnalysis(Atom atomHN, Atom atomN, Atom atomCA, Atom atomC, Atom atomO, Atom atomHA)
        {
            if ( atomN == null || atomN == null || atomCA == null || atomC == null || atomO == null || atomHA == null )
                throw new NullPointerException("missing atom in BackboneAnalysis");
            this.atomHN = atomHN;
            this.atomN = atomN;
            this.atomCA = atomCA;
            this.atomC = atomC;
            this.atomO = atomO;
            this.atomHA = atomHA;
        }
        /**
         * Static factory method to analyze the backbone atoms of a residue.
         * @param residue the residue to analyze
         */
        private static BackboneAnalysis getBackboneAnalysis(Peptide inputPeptide, Residue residue)
        {
            if ( ! inputPeptide.sequence.contains(residue) )
                throw new IllegalArgumentException("residue not in peptide");
            return new BackboneAnalysis(residue.HN, residue.N, residue.CA, residue.C, residue.O, residue.HA);
        }

        /**
         * Static factory method to analyze the backbone atoms of a ProtoAminoAcid.
         * @param targetPAA the ProtoAminoAcid to analyze
         */
        private static BackboneAnalysis getBackboneAnalysis(ProtoAminoAcid targetPAA)
        {
            Residue residue = targetPAA.residue;
            return new BackboneAnalysis(residue.HN, residue.N, residue.CA, residue.C, residue.O, residue.HA);
        }

        /**
         * Adjust the atom types of the input peitde to those in targetPAA.  The positions of the atoms will not
         * change.  The atom types are not adjusted if they are already the correct ones.  We clone the peptide with moveAtoms2
         * regardless of whether the atom types need adjusting because downstream operations will change the connectivity
         * graph of the peptide.
         *
         * Now also adjusts surface tensions through the Atom.changeTypes(atom) method.
         *
         * @return the adjusted peptide
         */
        public static Peptide adjustBackboneAtomTypes(Peptide inputPeptide, Residue inputResidue, ProtoAminoAcid targetPAA)
        {
            BackboneAnalysis residueBackbone = getBackboneAnalysis(inputPeptide, inputResidue);
            BackboneAnalysis protoAminoAcidBackbone = getBackboneAnalysis(targetPAA);
            
            // create a map from the old atoms to new atoms
            // the new atoms have the same positions but new atom types
            Map<Atom,Atom> atomMap = new HashMap<>();
            atomMap.put(residueBackbone.atomC, residueBackbone.atomC.changeTypes(protoAminoAcidBackbone.atomC));
            atomMap.put(residueBackbone.atomCA, residueBackbone.atomCA.changeTypes(protoAminoAcidBackbone.atomCA));
            atomMap.put(residueBackbone.atomN, residueBackbone.atomN.changeTypes(protoAminoAcidBackbone.atomN));
            atomMap.put(residueBackbone.atomO, residueBackbone.atomO.changeTypes(protoAminoAcidBackbone.atomO));
            
            // untested -- by analogy to HN
            if ( residueBackbone.atomHA != null && protoAminoAcidBackbone.atomHA != null )
                atomMap.put(residueBackbone.atomHA, residueBackbone.atomHA.changeTypes(protoAminoAcidBackbone.atomHA));
            
            // case 1: non-proline --> non-proline, HN is present in source and target, result: adjusted
            // case 2: non-proline -->     proline, HN is present in source but not in target,
            //                                      result: not adjusted because HN will be deleted
            // case 3: proline     --> non-proline, HN is not present in source but present in target,
            //                                      result: not adjusted because there is nothing to adjust yet
            // case 4: proline     --> proline,     HN is not present in source or target,
            //                                      result: not adjusted because there is nothing to adjust
            if (protoAminoAcidBackbone.atomHN != null && residueBackbone.atomHN != null)
                atomMap.put(residueBackbone.atomHN, residueBackbone.atomHN.changeTypes(protoAminoAcidBackbone.atomHN));
           
            // create the altered peptide
            return inputPeptide.moveAtoms2(atomMap);
        }

        @Override
        public int hashCode()
        {
            return Objects.hash(atomHN, atomN, atomCA, atomC, atomO, atomHA);
        }

        @Override
        public boolean equals(Object obj)
        {
            if ( obj == null )
                return false;
            if ( obj == this )
                return true;
            if ( !(obj instanceof BackboneAnalysis) )
                return false;

            BackboneAnalysis b = (BackboneAnalysis)obj;
            if ( Objects.equals(atomHN, b.atomHN) &&
                 Objects.equals(atomN,  b.atomN ) &&
                 Objects.equals(atomCA, b.atomCA) &&
                 Objects.equals(atomC,  b.atomC ) &&
                 Objects.equals(atomO,  b.atomO ) &&
                 Objects.equals(atomHA, b.atomHA)    )
                return true;
            return false;
        }
    }

    /**
     * Represents a new sidechain.
     */
    public static class Sidechain extends Molecule
    {
        /** for serialization */
        public static final long serialVersionUID = 1L;

        /** the beta carbon -- this atom will be connected to the backbone CA */
        public final Atom atomCB;

        /** the delta carbon if this is a proline -- this atom will be connected to the backbone N */
        public final Atom atomCD;

        /** the list of atoms that make up the chi angles in a side chain (i.e. CB, CG, CD, etc.) */
        public final List<Atom> chiAtoms;

        /**
         * This private constructor just copies fields.
         * @param atomCB the atom that should be attached to the backbone CA (cannot be null)
         * @param atomCD if non-null, the atom that should be attached to 
         * @param chiAtoms the atoms in the sidechain used for chi angle calculations
         */
        private Sidechain(String name, List<Atom> contents, SimpleWeightedGraph<Atom,DefaultWeightedEdge> connectivity,
                          Atom atomCB, Atom atomCD, List<Atom> chiAtoms)
        {
            super(name, contents, connectivity);
            if ( atomCB == null )
                throw new NullPointerException("null connection to backbone");
            this.atomCB = atomCB;
            // this should only be here in the case of a proline
            this.atomCD = atomCD;
            this.chiAtoms = chiAtoms;
        }

        /**
         * Makes the new sidechain for the mutated peptide.  Takes the sidechain template and moves it so that its CA-CB
         * vector will point in the same direction as the original CA-CB vector.  This first involves translation of the
         * sidechain atoms so that they are centered on CA.  Then, if necessary, a rotation to match the CA-CB vectors is
         * applied.  The new CA-CB distance is taken as the distance in the targetPAA template.
         * @return the new sidechain
         */
        @SuppressWarnings("unchecked")
        public static Sidechain makeNewSidechain(Peptide peptide, int residueIndex, ProtoAminoAcid targetPAA)
        {
            // get fields
            Residue inputResidue = peptide.sequence.get(residueIndex);

            // get the chi atoms before the sidechain has been moved to its new location
            List<Atom> beforeChiAtoms = new ArrayList<>();
            List<ProtoTorsion> templateChis = targetPAA.residue.chis;
            if ( templateChis.size() >= 1 )
                {
                    ProtoTorsion chi1 = templateChis.get(0);
                    beforeChiAtoms.add(chi1.atom3);
                    beforeChiAtoms.add(chi1.atom4);
                }
            if ( templateChis.size() >= 2 )
                {
                    for (int i=1; i < templateChis.size(); i++)
                        {
                            ProtoTorsion chiN = templateChis.get(i);
                            beforeChiAtoms.add(chiN.atom4);
                        }
                }
        
            // special fix for proline: deletes the N here because it's not part of the sidechain
            if ( targetPAA.residue.aminoAcid.isProline() )
                beforeChiAtoms.remove(beforeChiAtoms.size()-1);

            // get the sidechain atoms of the desired target
            // first, clone molecule
            Molecule targetPAAmolecule = targetPAA.molecule.moveAtoms(new HashMap<Atom,Atom>()); 
            
            // remove the ring connection if this is a proline
            Atom originalCD = null;
            if (targetPAA.residue.aminoAcid.isProline())
            {
                ProtoTorsion chi3 = templateChis.get(2);
                Atom atom3 = chi3.atom3;
                Atom atom4 = chi3.atom4;
                
                SimpleWeightedGraph<Atom,DefaultWeightedEdge> tempConnectivity = (SimpleWeightedGraph<Atom,DefaultWeightedEdge>)targetPAA.molecule.connectivity.clone();
                DefaultWeightedEdge e = tempConnectivity.removeEdge(atom3,atom4);
                if ( e == null)
                    throw new NullPointerException("unexpected null edge for proline connection");
                targetPAAmolecule = new Molecule(targetPAA.molecule.name, targetPAA.molecule.contents, tempConnectivity);
            
                for (Atom a : targetPAAmolecule.contents)
                    {
                        if ( a.type1 == 59 )
                            {
                                if ( originalCD == null )
                                    originalCD = a;
                                else
                                    throw new IllegalArgumentException("CD already set");
                            }
                    }
                if ( originalCD == null )
                    throw new NullPointerException("CD not found in a proline");
            }

            // get CA and CB in the source and target molecules
            Atom sourceCA = inputResidue.prochiralConnection.getFirst();
            Atom sourceCB = inputResidue.prochiralConnection.getSecond();
            Atom targetCA = targetPAA.residue.prochiralConnection.getFirst();
            Atom targetCB = targetPAA.residue.prochiralConnection.getSecond();
            if ( sourceCA == null || sourceCB == null || targetCA == null || targetCB == null )
                throw new NullPointerException("missing CA or CB");
            
            // translate the template sidechain so its CA is at the origin
            Vector3D targetTranslation = targetCA.position.negate();

            List<Atom> sidechainAtoms = new ArrayList<>(targetPAAmolecule.getHalfGraph(targetCA, targetCB));
            Map<Atom,Atom> atomMap1 = new HashMap<>();  // maps old template atoms to new template atoms
            
            for (Atom a : sidechainAtoms)
            {
                Vector3D newPosition = a.position.add(targetTranslation);
                Atom movedAtom = a.moveAtom(newPosition);
                atomMap1.put(a,movedAtom);
            }
           
            Vector3D newTemplateCBposition = targetCB.position.add(targetTranslation);

            // find the rotation that will make the new CA,CB vector match the original CA,CB vector
            // check for a degenerate rotation
            Vector3D sourceTranslation = sourceCA.position.negate();
            Vector3D translatedSourceCBposition = sourceCB.position.add(sourceTranslation);
            
            Vector3D rotationAxis = Vector3D.crossProduct(newTemplateCBposition, translatedSourceCBposition);
            Rotation rotation = null;
            try { rotation = new Rotation(rotationAxis, Vector3D.angle(translatedSourceCBposition, newTemplateCBposition)); }
            catch (Exception e) { rotation = Rotation.IDENTITY; }

            // apply the rotation and undo the translation
            Map<Atom,Atom> atomMap2 = new HashMap<>();
            for (Atom a : atomMap1.keySet())
                {
                    Atom oldAtom = atomMap1.get(a);
                    Vector3D oldPosition = oldAtom.position;
                    Vector3D newPosition = rotation.applyTo(oldPosition);
                    newPosition = newPosition.add(sourceTranslation.negate());
                    Atom rotatedAtom = oldAtom.moveAtom(newPosition);
                    atomMap2.put(oldAtom,rotatedAtom);
                }

            // update the sidechain atoms
            for (int i=0; i < sidechainAtoms.size(); i++)
                {
                    Atom oldAtom = sidechainAtoms.get(i);
                    Atom newAtom = atomMap1.get(oldAtom);
                    newAtom      = atomMap2.get(newAtom);
                    if ( newAtom == null )
                        throw new NullPointerException("null atom");
                    sidechainAtoms.set(i,newAtom);
                }

            // create a new connectivity graph
            SimpleWeightedGraph<Atom,DefaultWeightedEdge> sidechainConnectivity = new SimpleWeightedGraph<Atom,DefaultWeightedEdge>(DefaultWeightedEdge.class);
            for (Atom a : sidechainAtoms)
                sidechainConnectivity.addVertex(a);

            // add the old bonds
            for (DefaultWeightedEdge e : targetPAAmolecule.connectivity.edgeSet())
                {
                    Atom oldEdgeSource = targetPAAmolecule.connectivity.getEdgeSource(e);
                    Atom oldEdgeTarget = targetPAAmolecule.connectivity.getEdgeTarget(e);

                    Atom newEdgeSource = atomMap1.get(oldEdgeSource);
                    if ( newEdgeSource != null )
                        newEdgeSource = atomMap2.get(newEdgeSource);

                    Atom newEdgeTarget = atomMap1.get(oldEdgeTarget);
                    if ( newEdgeTarget != null )
                        newEdgeTarget = atomMap2.get(newEdgeTarget);

                    if ( newEdgeSource != null && newEdgeTarget != null )
                        sidechainConnectivity.addEdge(newEdgeSource, newEdgeTarget);
                }

            // locate the place where the sidechain will attach to the backbone (the beta carbon)
            Atom newCB = atomMap2.get(atomMap1.get(targetCB));
            Atom newCD = null;
            if ( originalCD != null )
                {
                    newCD = atomMap2.get(atomMap1.get(originalCD));
                    if ( newCD == null )
                        throw new NullPointerException("failure to set proline CD connection");
                }

            // update the chi atoms
            List<Atom> afterChiAtoms = new ArrayList<>();
            for (Atom oldAtom : beforeChiAtoms)
                {
                    Atom newAtom = atomMap1.get(oldAtom);
                    newAtom = atomMap2.get(newAtom);
                    if ( newAtom == null )
                        throw new NullPointerException("atom not found");
                    afterChiAtoms.add(newAtom);
                }
            
            // return the new Sidechain object
            return new Sidechain("sidechain", sidechainAtoms, sidechainConnectivity, newCB, newCD, afterChiAtoms);
        }
    }

    /** for testing */
    public static void main(String[] args)
    {
        DatabaseLoader.go();
        List<ProtoAminoAcid> sequence = ProtoAminoAcidDatabase.getSpecificSequence("l_pro","met","standard_ala","gly","d_proline", "gly", "phe", "val", "hd",         "l_pro");
        Peptide peptide = PeptideFactory.createPeptide(sequence);

        for (int i=0; i < peptide.sequence.size(); i++)
            {
                peptide = BackboneMutator.mutateOmega(peptide, i);
                peptide = BackboneMutator.mutatePhiPsi(peptide, i);
                peptide = RotamerMutator.mutateChis(peptide, i);
            }
        peptide = PeptideFactory.setHairpinAngles(peptide);
        System.out.println(peptide.sequence.get(0).toString(peptide));
        for (Atom a : peptide.sequence.get(0).atoms)
            System.out.println(a.toFullString());

        GaussianInputFile f = new GaussianInputFile(peptide);
        f.write("test_peptides/test.gjf");

        ProtoAminoAcid templatePAA = ProtoAminoAcidDatabase.getTemplate("standard_leucine");
        peptide = SidechainMutator.mutateSidechain(peptide, peptide.sequence.get(0), templatePAA);
        System.out.println(peptide.sequence.get(0).toString(peptide));
        for (Atom a : peptide.sequence.get(0).atoms)
            System.out.println(a.toFullString());

        f = new GaussianInputFile(peptide);
        f.write("test_peptides/test2.gjf");
    }
}
