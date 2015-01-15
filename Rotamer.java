import java.util.*;
import java.io.*;
import org.apache.commons.math3.geometry.euclidean.threed.*;
import org.jgrapht.*;
import org.jgrapht.graph.*;
import org.jgrapht.alg.*;
import com.google.common.collect.*;
import java.util.concurrent.*;

/**
 * Represents a sidechain at a particular position in a given peptide.
 * Note: if we start to do amino acids with two substituents instead of one, we'll run into a problem with
 * the way rotamer is defined, since a mutation from one amino acid to another won't preserve the backbone.
 */
public class Rotamer implements Immutable
{
    /** the atoms in the sidechain */
    public final List<Atom> atoms;

    /** the position in the sequence of the original peptide */
    public final int sequenceIndex;

    /** the chi angles this rotamer has */
    public final List<Double> chis;

    /** the ProtoAminoAcid template for this rotamer */
    public final ProtoAminoAcid protoAminoAcid;

    /**
     * records whether each atom is next to a polar atom
     * index by atom
     * the adjacent polar atom if it is
     * null if it isn't next to any polar atoms
     */
    public final List<Atom> nextToPolar;
    
    /** records the number of bonded neighbors for each atom (note won't deal with CA properly, but this shouldn't matter) */
    public final List<Integer> numberOfNeighbors;

    /**
     * The minimum probability of a rotamer before it will be considered.
     * {@link Settings#ROTAMER_LIBRARY_THRESHOLD} applies first.
     */
    public static final double PROBABILITY_THRESHOLD = 0.01;

    /** minimum threshold area for a peak to be used */
    public static final double PEAK_THRESHOLD = 0.025;

    /** the area of a peak is considered the integral over PEAK_SIZE points before and after the maximum */
    public static final int PEAK_SIZE = 3;

    /** the maximum number of peaks that will be returned */
    public static final int MAX_PEAKS = 2;

    /**
     * The maximum number of normal rotamers to return per position.  For normal rotameric amino acids,
     * this is what is sounds like.  For non-rotameric amino acids, we could get up to MAX_PEAKS *
     * MAX_ROTAMERS_PER_POSITION rotamers.
     */
    public static final int MAX_ROTAMERS_PER_POSITION = 10;

    /** elements that might serve as X in a hydrogen bond */
    public static final List<Element> HBOND_ELEMENTS = ImmutableList.of(Element.NITROGEN, Element.OXYGEN, Element.SULFUR);

    /** minimum X...H-X angle for a hydrogen bond in degrees */
    public static final double HBONDED_MINIMUM_ANGLE = 120.0;

    /** static initializer */
    static
    {
            // check class invariants
            if ( PEAK_SIZE < 0 || MAX_PEAKS < 1 )
                throw new IllegalArgumentException("check MAX_PEAKS / PEAK_SIZE");

            if ( PEAK_THRESHOLD < 0.0 || PEAK_THRESHOLD > 1.0 )
                throw new IllegalArgumentException("check PEAK_THRESHOLD");

            if ( PROBABILITY_THRESHOLD < 0.0 || PROBABILITY_THRESHOLD > 1.0 )    
                throw new IllegalArgumentException("check PROBABILITY_THRESHOLD");
    }

    /**
     * Constructor.  Just copies the fields in the constructor.
     */
    public Rotamer(List<Atom> atoms, int sequenceIndex, List<Double> chis, ProtoAminoAcid protoAminoAcid,
                    List<Atom> nextToPolar, List<Integer> numberOfNeighbors)
    {
        this.atoms = ImmutableList.copyOf(atoms);
        this.sequenceIndex = sequenceIndex;
        this.chis = ImmutableList.copyOf(chis);
        this.protoAminoAcid = protoAminoAcid;
        
        // necessary because ImmutableLists don't take nulls
        this.nextToPolar = Collections.unmodifiableList(nextToPolar);

        this.numberOfNeighbors = ImmutableList.copyOf(numberOfNeighbors);
    }

    /**
     * Finds the possible rotamers for a particular non-hairpin position of the peptide sequence.
     * Assumes that the residue is already of the correct type.
     * @param peptide the peptide to analyze
     * @param i the sequence index to analyze
     * @param backboneNHs the backbone NHs: Ns mapped to Hs
     * @param rotamerSpace the possible rotamers at every position (expected to be null when setting the transition states, but needed for histidines)
     * @param backboneAtoms the backbone atoms
     * @param possibleChis the possible chi angles for this residue (if set to null, chis will be automatically generated)
     * @return a list of rotamers that are possible at the specified position
     */
    @SuppressWarnings("unchecked")
    public static List<Rotamer> generateRotamers(Peptide peptide, int i, Set<Pair<Atom,Atom>> backboneNHs,
                                                 List<List<Rotamer>> rotamerSpace, List<Atom> backboneAtoms, List<List<Double>> possibleChis)
    {
        //if ( peptide.sequence.get(i).protoAminoAcid.r.aminoAcid.isProline() )
        //    System.out.println(" @ " + peptide.sequence.get(i).protoAminoAcid.hashCode());

        // this is the list we will be returning
        List<Rotamer> returnList = new ArrayList<>();

        // get the possible chi angles if necessary
        Residue r = peptide.sequence.get(i);
        AminoAcid aminoAcid = r.aminoAcid;
        AminoAcid.RotamerType rotamerType = r.aminoAcid.rotamerType;
        List<List<Double>> chis = null;
        if ( possibleChis != null )
            chis = possibleChis;
        else
            {
                if ( r.isHairpin )
                    {
                        // this residue should not be moved so return an empty list
                        //System.out.println("[]");
                        return returnList;
                    }
                else if ( rotamerType == AminoAcid.RotamerType.SPECIAL )
                    throw new IllegalArgumentException("special amino acid?");
                else
                    {
                        // fill a list with lists of possible chi angles
                        chis = getPossibleRotamers(peptide, r, backboneNHs, rotamerSpace, backboneAtoms);
                        if ( chis.size() == 0 && r.chis.size() > 0 )
                            return returnList;
                    }
            }

        // print debug information
        //System.out.println(r.description);
        //for (int j=0; j < chis.size(); j++)
        //    System.out.println("Rotamer " + j + ": " + chis.get(j).toString());

        // get the torsions that describe where the sidechain is
        LinkedList<ProtoTorsion> oldProtoTorsions = new LinkedList<>(r.chis);

        // get the connectivity graph
        // if this is a proline, make a copy of the graph and disconnect the chi3 bond to the backbone
        SimpleWeightedGraph<Atom,DefaultWeightedEdge> connectivity = peptide.connectivity;
        if ( aminoAcid.isProline() )
            {
                // clone the graph
                connectivity = (SimpleWeightedGraph<Atom,DefaultWeightedEdge>)connectivity.clone();

                // remove the bond that connects the proline ring to the backbone at the beta carbon
                ProtoTorsion chi3 = r.chis.get(2);
                Atom atom3 = chi3.atom3;
                Atom atom4 = chi3.atom4;
                DefaultWeightedEdge e = connectivity.removeEdge(atom3,atom4);
                if ( e == null )
                    throw new NullPointerException("unexpected null edge for proline connection");
            
                // don't set the last chi angle
                oldProtoTorsions.removeLast();
            }
        
        // get the atoms in the sidechain
        Pair<Atom,Atom> prochiralConnection = r.prochiralConnection;
        List<Atom> allAtoms = new ArrayList<>(getHalfGraph(connectivity, prochiralConnection.getFirst(),prochiralConnection.getSecond()));
    
        // populate fields that will be common to all the rotamer objects
        int sequenceIndex = i;
        ProtoAminoAcid protoAminoAcid = r.protoAminoAcid;
        //if ( protoAminoAcid.r.aminoAcid.isProline() )
        //    System.out.println("> " + protoAminoAcid.r.description + " : " + protoAminoAcid.hashCode());
        
        // for each atom, get the indices of all of its adjacent atoms and put them in a list
        // note: proline will not be handled correctly, but hydrogen bonds are not expected to it
        List<List<Integer>> adjacentAtomIndices = new ArrayList<>();
        List<Integer> numberOfNeighbors = new ArrayList<>();
        for (Atom a : allAtoms)
            {
                Set<Atom> adjacentAtoms = getAdjacentAtoms(peptide.connectivity,a);
                numberOfNeighbors.add(adjacentAtoms.size());
                List<Integer> list = new ArrayList<>();
                for (Atom a2 : adjacentAtoms)
                    {
                        int index = allAtoms.indexOf(a2);
                        // could return atoms outside of the sidechain, but only
                        // return the ones inside the sidechain
                        if ( index != -1 )
                            list.add(index);
                    }
                adjacentAtomIndices.add(list);
            }

        // add the backbone atoms we need for the first torsions
        if ( rotamerType == AminoAcid.RotamerType.HAS_NO_ROTAMERS )
            {
                // special case for where the amino acid has no rotamers
                List<Double> requestedAngles = ImmutableList.of();  // no torsion angles
                List<Atom> nextToPolar = new LinkedList<>();        // glycine and alanine don't have any polar atoms in their sidechains
                for (Atom a : allAtoms)
                    nextToPolar.add(null);
                //System.out.println(protoAminoAcid.r.description);
                Rotamer singleRotamer = new Rotamer(allAtoms, sequenceIndex, requestedAngles,
                                                    protoAminoAcid, nextToPolar, numberOfNeighbors);
                returnList.add(singleRotamer);
                return returnList;
            }

        allAtoms.add(oldProtoTorsions.get(0).atom1);
        allAtoms.add(oldProtoTorsions.get(0).atom2);

        // cycle through all the desired torsion angles, not moving anything, but collecting the indices of what to move
        Map<ProtoTorsion,List<Integer>> atomsToMoveMap = new HashMap<>();
        Map<ProtoTorsion,Integer> atom1indices = new HashMap<>();
        Map<ProtoTorsion,Integer> atom2indices = new HashMap<>();
        Map<ProtoTorsion,Integer> atom3indices = new HashMap<>();
        Map<ProtoTorsion,Integer> atom4indices = new HashMap<>();
        for (ProtoTorsion torsion : oldProtoTorsions)
            {
                // make a list of the atom indices to move
                // indices are defined as the indices (0...N-1) in allAtoms
                List<Integer> atomIndicesToMove = new LinkedList<>();
                
                // get the atoms to move from the connectivity graph
                // and add the indices (as defined in allAtoms) to a list
                // this way, we can look up a prototorsion and then see what atom indices we need to move
                //Set<Atom> check = peptide.getHalfGraph(torsion.atom3, torsion.atom4);
                //System.out.println(check.size());
                Set<Atom> atomsToMove = getHalfGraph(connectivity, torsion.atom2, torsion.atom3);

                //System.out.println(atomsToMove.size());
                for (Atom a : atomsToMove)
                    atomIndicesToMove.add(allAtoms.indexOf(a));

                atomsToMoveMap.put(torsion,atomIndicesToMove);

                // do the same thing for the prototorsion atoms so after each time we move a torsion
                // we can reconsitute the prototorsion
                int atom1index = allAtoms.indexOf(torsion.atom1);
                int atom2index = allAtoms.indexOf(torsion.atom2);
                int atom3index = allAtoms.indexOf(torsion.atom3);
                int atom4index = allAtoms.indexOf(torsion.atom4);

                // must deal with atoms that are part of the torsion but won't be in allAtoms

                atom1indices.put(torsion, atom1index);
                atom2indices.put(torsion, atom2index);
                atom3indices.put(torsion, atom3index);
                atom4indices.put(torsion, atom4index);
            }

        // cycle through all the desired torsion angles, actually moving the atoms this time
        for (int j=0; j < chis.size(); j++)
            {
                List<Double> requestedAngles = chis.get(j);

                // make a new list of atoms
                List<Atom> newAtoms = new ArrayList<>(allAtoms);

                // check the number of torsions we have is the same as the 
                //if ( oldProtoTorsions.size() != requestedAngles.size() )
                //    throw new IllegalArgumentException("torsion angle list size mismatch");

                // set the chi angles one after the other
                for (int k=0; k < oldProtoTorsions.size(); k++)
                    {
                        // get the old prototorsion
                        ProtoTorsion oldProtoTorsion = oldProtoTorsions.get(k);

                        // make a new prototorsion
                        // this is necessary because the atoms are changing between successive operations
                        Atom atom1 = newAtoms.get(atom1indices.get(oldProtoTorsion));
                        Atom atom2 = newAtoms.get(atom2indices.get(oldProtoTorsion));
                        Atom atom3 = newAtoms.get(atom3indices.get(oldProtoTorsion));
                        Atom atom4 = newAtoms.get(atom4indices.get(oldProtoTorsion));
                        ProtoTorsion newProtoTorsion = new ProtoTorsion(atom1, atom2, atom3, atom4);

                        // figure out what angle we want to set the dihedral to
                        double theta = requestedAngles.get(k);

                        // figure out which atoms we have to move
                        Set<Atom> atomsToMove = new HashSet<>();
                        for (Integer l : atomsToMoveMap.get(oldProtoTorsion))
                            atomsToMove.add(newAtoms.get(l));

                        // move the atoms
                        newAtoms = setDihedral(newAtoms, atomsToMove, newProtoTorsion, theta);
                    }

                // remove the extra atoms
                newAtoms.remove(newAtoms.size()-1);
                newAtoms.remove(newAtoms.size()-1);

                // add the adjacent polar atoms
                List<Atom> nextToPolar = new ArrayList<>();  // parallel list to newAtoms, the adjacent polar atom if present, null otherwise
                for (int k=0; k < newAtoms.size(); k++)
                    {
                        Atom thisAtom = newAtoms.get(k);

                        // only note any adjacent polar atoms if this is a hydrogen atom
                        if ( thisAtom.element != Element.HYDROGEN )
                            {
                                nextToPolar.add(null);
                                continue;
                            }

                        // get the indices and translate back into atoms
                        List<Integer> list = adjacentAtomIndices.get(k);
                        Atom polarAtom = null;
                        for (Integer index : list)
                            {
                                Atom a = newAtoms.get(index);
                                if ( HBOND_ELEMENTS.contains(a.element) )
                                    {
                                        polarAtom = a;
                                        break;
                                    }
                            }
                        nextToPolar.add(polarAtom);
                    }

                // check
                //int expectedSize = newAtoms.size();
                //if ( nextToPolar.size() != expectedSize || numberOfNeighbors.size() != expectedSize )
                //    throw new IllegalArgumentException("list size mismatch");

                // return the new Rotamer
                Rotamer newRotamer = new Rotamer(newAtoms, sequenceIndex, requestedAngles,
                                                 protoAminoAcid, nextToPolar, numberOfNeighbors);
                returnList.add(newRotamer);
            }

        // return the result
        return returnList;
    }

    /**
     * Mutates the given peptide to the specified rotamer. This will automatically
     * perform the mutation at the correct location in the original peptide. If the identity of the residue in
     * startingPeptide is wrong, we fix it first. We don't check that the rotamer is appropriate for this position.
     * @param startingPeptide the peptide to make the mutation on
     * @param rotamer the rotamer we want in the resulting peptide
     * @return the mutated peptide
     */
    public static Peptide reconstitute(Peptide startingPeptide, Rotamer rotamer)
    {
        // check if the starting peptide requires a sidechain mutation
        Residue originalResidue = startingPeptide.sequence.get(rotamer.sequenceIndex);
        Peptide currentPeptide = startingPeptide;
        if ( ! originalResidue.protoAminoAcid.equals(rotamer.protoAminoAcid) )
            currentPeptide = SidechainMutator.mutateSidechain(currentPeptide, originalResidue, rotamer.protoAminoAcid);

        // adjust the chis
        if (rotamer.chis.size() > 0)
            currentPeptide = RotamerMutator.setChis(currentPeptide, currentPeptide.sequence.get(rotamer.sequenceIndex), rotamer.chis);
        return currentPeptide;
    }

    /**
     * Mutates the given peptide to the specified rotamers.  No checks.
     * @param startingPeptide the starting peptide
     * @param rotamerList the rotamers to use
     * @return the mutated peptide
     */
    public static Peptide reconstitute(Peptide startingPeptide, List<Rotamer> rotamerList)
    {
        Peptide returnPeptide = startingPeptide;
        for (Rotamer r : rotamerList)
            returnPeptide = reconstitute(returnPeptide, r);
        return returnPeptide;
    }

    /**
     * Rotates the atoms in the specified list by the desired amount.
     * The advantage in using this over the related method {@link Molecule#setDihedral(ProtoTorsion,double)} is that
     * it's faster.  We don't have to rebuild the entire connectivity graph.  We need a ProtoTorsion so we know
     * how much to rotate by, but we don't have to search the entire connectivity graph over and over
     * for the atoms we need to rotate.  The order of the atoms in the original list is preserved.
     * @param allAtoms all of the atoms in this fragment
     * @param atomsToRotate the atoms we should be rotating
     * @param protoTorsion the torsion to rotate
     * @param theta the angle we want to end up with (in degrees)
     * @return the rotated fragment
     */
    public static List<Atom> setDihedral(List<Atom> allAtoms, Collection<Atom> atomsToRotate, ProtoTorsion protoTorsion, double theta)
    {
        // determine how much rotation is needed
        double currentDihedralAngle = protoTorsion.getDihedralAngle();
        double requiredRotation = currentDihedralAngle - theta;

        // get the indices of the atoms to rotate
        List<Integer> indicesToRotate = new LinkedList<>();
        for (Atom a : atomsToRotate)
            {
                int index = allAtoms.indexOf(a);
                if ( index == -1 )
                    throw new IllegalArgumentException("atom to rotate is not in allAtoms");
                indicesToRotate.add(index);
            }

        // determine the rotation axis and create the rotation
        Vector3D atom2position = protoTorsion.atom2.position;
        Vector3D atom3position = protoTorsion.atom3.position;
        Vector3D rotationAxis = atom2position.subtract(atom3position);
        Rotation rotation = new Rotation(rotationAxis, Math.toRadians(requiredRotation));

        // perform the rotation
        List<Atom> returnList = new ArrayList<>(allAtoms.size());
        for (int i=0; i < allAtoms.size(); i++)
            {
                Atom a = allAtoms.get(i);
                if ( indicesToRotate.contains(i) )
                    {
                        // set the origin to atom3 by translating
                        Vector3D newPosition = a.position.subtract(atom3position);

                        // apply the rotation
                        newPosition = rotation.applyTo(newPosition);

                        // undo the translation
                        newPosition = newPosition.add(atom3position);

                        // apply the result
                        Atom newAtom = a.moveAtom(newPosition);
                        returnList.add(newAtom);
                    }
                else
                    returnList.add(a);
            }

        // return the result
        return returnList;
    }

    /**
     * Given a bond between includeAtom and excludeAtom, returns a list of atoms
     * containing all the atoms on the includeAtom side of the bond, including
     * includeAtom.  Will throw an exception if includeAtom and excludeAtom
     * form a ring or either atom is not in the connectivity graph.
     * @param connectivity the connectivity graph to search
     * @param excludeAtom this atom will not be included in the result
     * @param includeAtom this atom will be included in the result
     * @return the atoms on the includeAtom side of the graph
     */
    public static Set<Atom> getHalfGraph(SimpleWeightedGraph<Atom,DefaultWeightedEdge> connectivity, Atom excludeAtom, Atom includeAtom)
    {
        Set<Atom> returnSet = new HashSet<Atom>();

        // check conditions
        if ( ! connectivity.containsVertex(excludeAtom) )
            throw new IllegalArgumentException("exclude atom not in graph");
        if ( ! connectivity.containsVertex(includeAtom) )
            throw new IllegalArgumentException("include atom not in graph");
        
        if ( ! getAdjacentAtoms(connectivity,excludeAtom).contains(includeAtom) )
            throw new IllegalArgumentException("not bonded");

        // preform a breadth-first search of one branch of the graph only
        LinkedList<Atom> searchQueue = new LinkedList<Atom>();

        for (Atom a : getAdjacentAtoms(connectivity,includeAtom))
            {
                searchQueue.add(a);
                returnSet.add(a);
            }
        searchQueue.remove(includeAtom);
        searchQueue.remove(excludeAtom);
        returnSet.remove(includeAtom);
        returnSet.remove(excludeAtom);

        while (searchQueue.size() > 0)
            {
                Atom currentNode = searchQueue.remove();
                
                for (Atom a : getAdjacentAtoms(connectivity,currentNode))
                    {
                        // if the excluded atom is found, this is a ring!
                        if ( a == excludeAtom )
                            {
                                System.out.println(excludeAtom);
                                System.out.println(includeAtom);
                                throw new IllegalArgumentException("found a ring during half graph search");
                            }

                        // if a isn't in returnSet, add it, and queue it for investigation
                        if ( ! returnSet.contains(a) && a != includeAtom )
                            {
                                returnSet.add(a);
                                searchQueue.add(a);
                            }
                    }
            }
        returnSet.add(includeAtom);
        return(returnSet);
    }

    /**
     * Returns the directly bonded neighbors of the specified atom.
     * @param connectivity the connectivity graph to search
     * @param atom the target vertex to search
     */
    public static Set<Atom> getAdjacentAtoms(SimpleWeightedGraph<Atom,DefaultWeightedEdge> connectivity, Atom atom)
    {
        Set<Atom> returnSet = new HashSet<Atom>();
        if ( ! connectivity.containsVertex(atom) )
            throw new IllegalArgumentException("atom must be within this connectivity graph!");
        for (DefaultWeightedEdge e : connectivity.edgesOf(atom))
            {
                returnSet.add(connectivity.getEdgeSource(e));
                returnSet.add(connectivity.getEdgeTarget(e));
            }
        returnSet.remove(atom);
        return(returnSet);
    }

    /** how many points to use to search for TS rotamers */
    public static final int TS_GRID_SIZE = 13;

    /**
     * Finds the chi values for a transition state amino acid, under the constraint that the
     * anionic oxygen must be near a backbone HN.  Clashes with the backbone are crudely checked.
     * @param peptide the peptide
     * @param residue the residue where the transition state currently is
     * @param backboneNHs backbone Ns mapped to backbone Hs
     * @param backboneAtoms the backbone atoms
     * @return lists of lists (chi1, chi2) for transition rotamers that satisfy the constraint
     */
    public static List<List<Double>> getTSchis(Peptide peptide, Residue residue, Set<Pair<Atom,Atom>> backboneNHs, List<Atom> backboneAtoms)
    {
        // get the rotation axes
        ProtoTorsion chi1torsion = residue.chis.get(0);
        ProtoTorsion chi2torsion = residue.chis.get(1);

        // get the sidechain positions
        Pair<Atom,Atom> prochiralConnection = residue.prochiralConnection;
        List<Atom> allAtoms = new ArrayList<>(peptide.getHalfGraph(prochiralConnection.getFirst(), prochiralConnection.getSecond()));
        allAtoms.add(chi1torsion.atom1);
        allAtoms.add(chi1torsion.atom2);

        // we don't want to check atoms that are part of the backbone for clashes
        // but just in case there are numerical issues we refer to them by index instead of reference
        HashSet<Integer> ignoreAtomIndices = new HashSet<>();
        ignoreAtomIndices.add(allAtoms.indexOf(chi1torsion.atom1));
        ignoreAtomIndices.add(allAtoms.indexOf(chi1torsion.atom2));
        ignoreAtomIndices.add(allAtoms.indexOf(chi1torsion.atom3));

        // the atoms to rotate for each torsion
        List<Atom> chi1atoms = new ArrayList<>(peptide.getHalfGraph(chi1torsion.atom2, chi1torsion.atom3));
        List<Atom> chi2atoms = new ArrayList<>(peptide.getHalfGraph(chi2torsion.atom2, chi2torsion.atom3));

        // the indices to rotate for chi2
        List<Integer> chi2indices = new ArrayList<>();
        for (Atom a : chi2atoms)
            {
                int atomIndex = allAtoms.indexOf(a);
                chi2indices.add(atomIndex);
            }

        // the indices of the torsion for chi2
        List<Integer> chi2torsionIndices = new ArrayList<>();
        chi2torsionIndices.add(allAtoms.indexOf(chi2torsion.atom1));
        chi2torsionIndices.add(allAtoms.indexOf(chi2torsion.atom2));
        chi2torsionIndices.add(allAtoms.indexOf(chi2torsion.atom3));
        chi2torsionIndices.add(allAtoms.indexOf(chi2torsion.atom4));

        // get the index of the transition state oxygen
        Integer atomOindex = null; 
        for (Atom a : allAtoms)
            {
                if ( a.tinkerAtomType == 408 )
                    {
                        atomOindex = allAtoms.indexOf(a);
                        break;
                    }
            }
        if ( atomOindex == null )
            throw new NullPointerException("null oxygen index");


        // rotate chi1 and chi2 on a grid
        // check for the desired contact
        List<List<Double>> returnList = new ArrayList<>();
        double stepSize = 360.0 / (TS_GRID_SIZE-1.0);
        double chi1 = -180.0;
        int sequenceIndex = peptide.sequence.indexOf(residue);
        if ( sequenceIndex == -1 )
            throw new IllegalArgumentException("bad sequence index");
        for (int i=0; i < TS_GRID_SIZE; i++)
            {
                List<Atom> thisAtoms = setDihedral(allAtoms, chi1atoms, chi1torsion, chi1);
                chi2atoms.clear();
                for (Integer index : chi2indices)
                    chi2atoms.add(thisAtoms.get(index));
                Atom atom1 = thisAtoms.get(chi2torsionIndices.get(0));
                Atom atom2 = thisAtoms.get(chi2torsionIndices.get(1));
                Atom atom3 = thisAtoms.get(chi2torsionIndices.get(2));
                Atom atom4 = thisAtoms.get(chi2torsionIndices.get(3));
                chi2torsion = new ProtoTorsion(atom1, atom2, atom3, atom4);

                double chi2 = -180.0;
                for (int j=0; j < TS_GRID_SIZE; j++)
                    {
                        // try this rotamer
                        List<Atom> thisAtoms2 = setDihedral(thisAtoms, chi2atoms, chi2torsion, chi2);

                        // get the oxygen
                        Atom atomO = thisAtoms2.get(atomOindex);

                        // check that the oxygen forms a hydrogen bond to a backbone NH
                        boolean acceptable = false;
                        for (Pair<Atom,Atom> pair : backboneNHs)
                            {
                                Atom atomN = pair.getFirst();
                                Atom atomH = pair.getSecond();

                                // check distance and angle
                                double distance = Vector3D.distance(atomO.position, atomH.position);
                                if ( distance > 2.50 )
                                    continue;
                                
                                double angle = Molecule.getAngle(atomO.position, atomH.position, atomN.position);
                                if ( distance > RotamerPacker.HBONDED_MIN_DISTANCE && angle > RotamerPacker.HBONDED_MINIMUM_ANGLE )
                                    {
                                        acceptable = true;
                                        break;
                                    }
                            }
                        
                        // if acceptable, add it to the list of chis to return
                        if ( acceptable )
                            {
                                // check for backbone atom clashes
                                HashSet<Atom> ignoreAtoms = new HashSet<>();
                                for (Integer index : ignoreAtomIndices)
                                    ignoreAtoms.add(thisAtoms2.get(index));

                                // check if this rotamer now clashes with the backbone
                                checking:
                                for (Atom rotamerAtom : thisAtoms2)
                                    {
                                        if ( ignoreAtoms.contains(rotamerAtom) )
                                            continue;

                                        for (Atom backboneAtom : backboneAtoms)
                                            {
                                                double distance = Vector3D.distance(rotamerAtom.position, backboneAtom.position);
                                                if ( distance < RotamerPacker.HBONDED_MIN_DISTANCE )
                                                    {
                                                        acceptable = false;
                                                        break checking;
                                                    }
                                            }
                                    }
                                if (acceptable)
                                    {
                                        List<Double> thisChis = ImmutableList.of(chi1,chi2);
                                        returnList.add(thisChis);
                                    }
                            }

                        // increment chi2
                        chi2 += stepSize;
                    }
                // increment chi1
                chi1 += stepSize;
            }
        //System.out.println("ts: " + returnList);
        return ImmutableList.copyOf(returnList);
    }

    /**
     * Finds the chi values for a histidine amino acid, under the constraint that the pi-nitrogen must be near a
     * transition state hydroxyl.
     */
    public static List<List<Double>> getHistidineChis(Peptide peptide, Residue residue,
                                                      List<List<Rotamer>> rotamerSpace, List<Atom> backboneAtoms)
    {
        // make a list of interesting atoms at other positions; i.e., the hydroxyl hydrogens
        int i = peptide.sequence.indexOf(residue);
        List<Atom> interestingAtoms = new ArrayList<>();
        List<Atom> interestingOxygenAtoms = new ArrayList<>(); // make a list of the corresponding hydroxyl oxygens
        for (int j=0; j < rotamerSpace.size(); j++)
            {
                if ( i==j )
                    continue;
                List<Rotamer> list = rotamerSpace.get(j);
                for (Rotamer r : list)
                    {
                        if ( r.protoAminoAcid.r.description.indexOf("transition_state") > -1 )
                            {
                                interestingAtoms.add(r.getInterestingAtom());
                                boolean partnerFound = false;
                                for (Atom a : r.atoms)
                                    {
                                        if ( a.tinkerAtomType == 402 )
                                            {
                                                interestingOxygenAtoms.add(a);
                                                partnerFound = true;
                                                break;
                                            }
                                    }
                                if ( !partnerFound )
                                    throw new IllegalArgumentException("partner not found");
                            }
                    }
            }


        // get the rotation axes
        ProtoTorsion chi1torsion = residue.chis.get(0);
        ProtoTorsion chi2torsion = residue.chis.get(1);

        // get the sidechain positions
        Pair<Atom,Atom> prochiralConnection = residue.prochiralConnection;
        List<Atom> allAtoms = new ArrayList<>(peptide.getHalfGraph(prochiralConnection.getFirst(), prochiralConnection.getSecond()));
        allAtoms.add(chi1torsion.atom1);
        allAtoms.add(chi1torsion.atom2);

        // we don't want to check atoms that are part of the backbone for clashes
        // but just in case there are numerical issues we refer to them by index instead of reference
        HashSet<Integer> ignoreAtomIndices = new HashSet<>();
        ignoreAtomIndices.add(allAtoms.indexOf(chi1torsion.atom1));
        ignoreAtomIndices.add(allAtoms.indexOf(chi1torsion.atom2));
        ignoreAtomIndices.add(allAtoms.indexOf(chi1torsion.atom3));
        
        // the atoms to rotate for each torsion
        List<Atom> chi1atoms = new ArrayList<>(peptide.getHalfGraph(chi1torsion.atom2, chi1torsion.atom3));
        List<Atom> chi2atoms = new ArrayList<>(peptide.getHalfGraph(chi2torsion.atom2, chi2torsion.atom3));

        // the indices to rotate for chi2
        List<Integer> chi2indices = new ArrayList<>();
        for (Atom a : chi2atoms)
            {
                int atomIndex = allAtoms.indexOf(a);
                chi2indices.add(atomIndex);
            }

        // the indices of the torsion for chi2
        List<Integer> chi2torsionIndices = new ArrayList<>();
        chi2torsionIndices.add(allAtoms.indexOf(chi2torsion.atom1));
        chi2torsionIndices.add(allAtoms.indexOf(chi2torsion.atom2));
        chi2torsionIndices.add(allAtoms.indexOf(chi2torsion.atom3));
        chi2torsionIndices.add(allAtoms.indexOf(chi2torsion.atom4));

        // get the index of the histidine pi-nitrogen
        Integer atomNindex = null; 
        for (Atom a : allAtoms)
            {
                if ( a.tinkerAtomType == 126 || a.tinkerAtomType == 135 )
                    {
                        atomNindex = allAtoms.indexOf(a);
                        break;
                    }
            }
        if ( atomNindex == null )
            throw new NullPointerException("null histidine N index");

        // rotate chi1 and chi2 on a grid
        // check for the desired contact
        List<List<Double>> returnList = new ArrayList<>();
        double stepSize = 360.0 / (TS_GRID_SIZE-1.0);
        double chi1 = -180.0;
        int sequenceIndex = peptide.sequence.indexOf(residue);
        if ( sequenceIndex == -1 )
            throw new IllegalArgumentException("bad sequence index");
        for (int j=0; j < TS_GRID_SIZE; j++)
            {
                List<Atom> thisAtoms = setDihedral(allAtoms, chi1atoms, chi1torsion, chi1);
                chi2atoms.clear();
                for (Integer index : chi2indices)
                    chi2atoms.add(thisAtoms.get(index));
                Atom atom1 = thisAtoms.get(chi2torsionIndices.get(0));
                Atom atom2 = thisAtoms.get(chi2torsionIndices.get(1));
                Atom atom3 = thisAtoms.get(chi2torsionIndices.get(2));
                Atom atom4 = thisAtoms.get(chi2torsionIndices.get(3));
                chi2torsion = new ProtoTorsion(atom1, atom2, atom3, atom4);

                double chi2 = -180.0;
                for (int k=0; k < TS_GRID_SIZE; k++)
                    {
                        // try this rotamer
                        List<Atom> thisAtoms2 = setDihedral(thisAtoms, chi2atoms, chi2torsion, chi2);

                        // get the oxygen
                        Atom atomN = thisAtoms2.get(atomNindex);

                        // check that this could become an interesting pair
                        boolean acceptable = false;

                        for ( Atom a : interestingAtoms )
                            {
                                double distance = Vector3D.distance(atomN.position, a.position);
                                if ( distance > RotamerPacker.HBONDED_MIN_DISTANCE && distance < RotamerPacker.INTERESTING_DISTANCE )
                                    {
                                        Integer index = interestingAtoms.indexOf(a); 
                                        Atom oxygenAtom = interestingOxygenAtoms.get(index);
                                        double angle = Molecule.getAngle(oxygenAtom, a, atomN); 
                                        if ( angle > HBONDED_MINIMUM_ANGLE )
                                            {
                                                acceptable = true;
                                                break;
                                            }
                                    }
                            }

                        if ( acceptable )
                            {
                                // check for backbone atom clashes
                                HashSet<Atom> ignoreAtoms = new HashSet<>();
                                for (Integer index : ignoreAtomIndices)
                                    ignoreAtoms.add(thisAtoms2.get(index));

                                clashCheck:
                                for ( Atom a1 : thisAtoms2 )
                                    {
                                        if ( ignoreAtoms.contains(a1) )
                                            continue;
                                        for (Atom a2 : backboneAtoms )
                                            {
                                                double distance = Vector3D.distance(a1.position, a2.position);
                                                if ( distance < RotamerPacker.HBONDED_MIN_DISTANCE )
                                                    {
                                                        acceptable = false;
                                                        break clashCheck;
                                                    }
                                            }
                                    }
                        
                                // if acceptable, add it to the list of chis to return
                                if ( acceptable )
                                    returnList.add(ImmutableList.of(chi1,chi2));
                            }

                        // increment chi2
                        chi2 += stepSize;
                    }
                // increment chi1
                chi1 += stepSize;
            }
        //System.out.println("his size: " + returnList.size());
        //System.out.println("his: " + returnList);
        return returnList;
    }

    /**
     * Given a residue, return the possible rotamers.  Rotamers are based on the current
     * phi and psi values.  These are only the rotamers for the current amino acid type.
     *
     * Will throw an exception if this residue should not be rotated.  Rotamers beneath
     * PROBABILITY_THRESHOLD will be ignored.  Non-rotameric degrees of freedom
     * will be summarized into several angles using BIN_SIZE and BIN_THRESHOLD.
     *
     * @param peptide the peptide the rotamer should be generated in
     * @param residue the input residue
     * @param backboneNHs the backbone NHs in this peptide
     * @param rotamerSpace the existing rotamer space
     * @param backboneAtoms the backbone atoms
     * @return a nested list of all the angles chi1, chi2, ..., chiN
     */
    public static List<List<Double>> getPossibleRotamers(Peptide peptide, Residue residue, Set<Pair<Atom,Atom>> backboneNHs,
                                                         List<List<Rotamer>> rotamerSpace, List<Atom> backboneAtoms)
    {
        AminoAcid.RotamerType rotamerType = residue.aminoAcid.rotamerType;
        AminoAcid aminoAcid = residue.aminoAcid;
        Map<AminoAcid,SidechainRotamerLibrary> MAP = RotamerLibrary.MAP;

        double phi = residue.phi.getDihedralAngle();
        double psi = residue.psi.getDihedralAngle();

        if ( residue.description.indexOf("hairpin") > -1 )
            throw new IllegalArgumentException("won't rotate hairpin residues");
        else if ( residue.description.indexOf("transition_state") > -1 )
            {
                // dynamically generate chis for transition state amino acids
                return getTSchis(peptide,residue,backboneNHs,backboneAtoms);
            }
        else if ( residue.description.indexOf("histidine") > -1 )
            {
                // dynamically generate chis for histidines based on what transition states are already present
                return getHistidineChis(peptide,residue,rotamerSpace, backboneAtoms);
            }
        else if (rotamerType == AminoAcid.RotamerType.HAS_NO_ROTAMERS)
            {
                //throw new IllegalArgumentException("non rotameric amino acids don't have rotamers (" + aminoAcid.toString() + ")");
                List<List<Double>> returnList = new LinkedList<>();
                returnList.add(new LinkedList<Double>());
                return returnList;
            }
        else if (rotamerType == AminoAcid.RotamerType.SPECIAL)
            throw new IllegalArgumentException("special amino acids don't have rotamers (" + aminoAcid.toString() + ")");
        else if (rotamerType == AminoAcid.RotamerType.IS_ROTAMERIC)
            {
                RotamericLibrary rotLib = (RotamericLibrary) MAP.get(aminoAcid);
                DiscreteProbabilityDistribution<List<Double>> dpd = rotLib.get(phi,psi);
                List<List<Double>> outcomes = dpd.outcomes;
                List<Double> probabilities = dpd.inputProbabilities;

                // sort in descending order
                TreeMap<Double,List<Double>> allRotamers = new TreeMap<>(Collections.reverseOrder());
                
                for (int i=0; i < outcomes.size(); i++)
                    {
                        List<Double> outcome = outcomes.get(i);
                        Double probability = probabilities.get(i);
                        if ( probability < PROBABILITY_THRESHOLD )
                            continue;

                        // map probability to rotamers so we can sort to get the highest probability rotamers
                        allRotamers.put(probability,outcome);
                    }
                if ( allRotamers.size() == 0 )
                    throw new IllegalArgumentException("expected to find rotamers (rotameric)");
                
                // only include up to MAX_ROTAMER_PER_POSITION rotamers in the final list
                List<List<Double>> prunedRotamers = new LinkedList<>();
                for (Double probability : allRotamers.keySet())
                    {
                        if ( prunedRotamers.size() > MAX_ROTAMERS_PER_POSITION )
                            break;
                        List<Double> thisRotamer = allRotamers.get(probability);
                        prunedRotamers.add(thisRotamer);
                    }
                
                return prunedRotamers;
            }

        else if (rotamerType == AminoAcid.RotamerType.NON_ROTAMERIC)
            {
                NonRotamericLibrary nRotLib = (NonRotamericLibrary) MAP.get(aminoAcid);
                DiscreteProbabilityDistribution<NonRotamericLibrary.NonRotamericAngles> dpd = nRotLib.get(phi,psi);
                List<NonRotamericLibrary.NonRotamericAngles> outcomes = dpd.outcomes;
                List<Double> probabilities = dpd.inputProbabilities;
                
                // find the most probable rotamers
                // create sorted map in descending order
                TreeMap<Double,NonRotamericLibrary.NonRotamericAngles> allRotamers = new TreeMap<>(Collections.reverseOrder());
                for (int i=0; i < outcomes.size(); i++)
                    {
                        NonRotamericLibrary.NonRotamericAngles outcome = outcomes.get(i);
                        Double probability = probabilities.get(i);
                        if ( probability < PROBABILITY_THRESHOLD )
                            continue;
                        allRotamers.put(probability,outcome);
                    }
                
                List<NonRotamericLibrary.NonRotamericAngles> prunedRotamers = new LinkedList<>();
                for (Double probability : allRotamers.keySet())
                    {
                        if (prunedRotamers.size() > MAX_ROTAMERS_PER_POSITION)
                            break;
                        NonRotamericLibrary.NonRotamericAngles outcome = allRotamers.get(probability);
                        prunedRotamers.add(outcome);
                    }

                // convert from non-rotameric angles to rotameric angles
                List<List<Double>> returnList = new LinkedList<>();
                for (NonRotamericLibrary.NonRotamericAngles outcome : prunedRotamers)
                    {
                        // convert outcome to rotamer angles in the form of List<Double>
                        // to do this, we have to deal with the non-rotameric degree of freedom
                        //
                        // first, get the standard rotameric angles and the 
                        // enclosed dpd for the standard chi angles
                        List<Double> rotamericAngles = outcome.getRotamericAngles();
                        DiscreteProbabilityDistribution<Double> dpd1 = outcome.getDPD();
                        
                        //System.out.println("phi: " + phi);
                        //System.out.println("psi: " + psi);
                        //System.out.println(rotamericAngles);

                        // get some representative values of the non-rotameric torsion angle
                        List<Double> nonRotamericAngles = summarize(dpd1);

                        // combine the rotameric angles with the summarized nonRotamericAngles
                        // to make several new overall rotamers
                        for (Double lastAngle : nonRotamericAngles)
                            {
                                List<Double> thisRotamer = new LinkedList<>(rotamericAngles);
                                thisRotamer.add(lastAngle);
                                returnList.add(thisRotamer);
                                //returnList.add(ImmutableList.copyOf(thisRotamer));
                            }
                    }

                if ( returnList.size() == 0 )
                    throw new IllegalArgumentException("expected to find rotamers (non-rotameric)");
                return returnList;
                //return ImmutableList.copyOf(returnList);
            }

        // should be unreachable
        return null;
    }

    /**
     * Finds the expected value of the given distribution (probability-weighted average
     * over all outcomes).  Only applies to distributions of doubles.
     * @param dpd the distribution
     * @return the expected value
     */
    public static Double getExpectedValue(DiscreteProbabilityDistribution<Double> dpd)
    {
        // get data
        List<Double> probabilities = new ArrayList<>(dpd.inputProbabilities);
        List<Double> outcomes = new ArrayList<>(dpd.outcomes);
        
        // discard duplicate data
        // (+180 degrees is the same as -180 degrees)
        boolean remove = false;
        for (Double d : outcomes)
            {
                if ( d == -180.0 )
                    remove = true;
            }

        if ( remove )
            {
                for (int i=0; i < probabilities.size(); i++)
                    {
                        Double outcome = outcomes.get(i);
                        if ( outcome == 180.0 )
                            {
                                probabilities.remove(i);
                                outcomes.remove(i);
                            }
                    }
            }

        // normalize
        probabilities = normalize(probabilities);

        // calculate expected value
        double expectedValue = 0.0;
        for (int i=0; i < probabilities.size(); i++)
            {
                Double probability = probabilities.get(i);
                Double outcome = outcomes.get(i);
                expectedValue += probability * outcome;
            }

        return expectedValue;
    }

    /**
     * Takes a list of non-negative doubles and normalizes it.  That is, each value will
     * be divided by the original sum of values to produce a list whose sum is 1.0.
     * An exception will be thrown if an input value is negative.  Returns a mutable
     * ArrayList.
     * @param list the input list containin non-negative doubles
     * @return the normalized list
     */
    public static ArrayList<Double> normalize(List<Double> list)
    {
        for ( Double d : list )
            {
                if ( d < 0.0 )
                    throw new IllegalArgumentException("negative numbers are not allowed");
            }
        ArrayList<Double> returnList = new ArrayList<>(list);
        
        // calculate sum of probabilities
        double sum = 0.0;
        for (Double d : returnList)
            sum += d;

        // normalize probabilities first, since the input probabilities might not be normalized
        for (int i=0; i < returnList.size(); i++)
            returnList.set(i, returnList.get(i) / sum);
        
        return returnList;
    }

    /**
     * Takes a histogram of double values and returns a small number of doubles
     * that are considered representative.  These values are selected by finding
     * the positions of peaks in the histogram.  The histogram is assumed to be
     * reasonably smooth and a peak is defined as a point whose neighbors are
     * lower in value.  Additionally, the area under the peak (PEAK_SIZE points to
     * the left and right) must exceed PEAK_THRESHOLD.  The peaks are sorted by
     * area and the peaks with the largest areas are considered.  (At most, MAX_PEAKS
     * peaks will be returned.)  The expected values of these peaks are calculated
     * and returned.
     *
     * Note -- this does not account for the fact that the interval -180,180 is cyclic.
     * Edge peaks might be counted twice, although I think it's unlikely.  Also, small
     * probability outcomes are pruned from the histogram coming in, so dpd might not
     * span the [-180,180] interval.
     *
     * @param dpd the input distribution
     * @return the summarized values
     */
    public static List<Double> summarize(DiscreteProbabilityDistribution<Double> dpd)
    {
        // obtain the underlying data
        List<Double> probabilities = normalize(dpd.inputProbabilities);
        List<Double> outcomes = dpd.outcomes;
        //System.out.println("expected value: " + getExpectedValue(dpd));
        //writeCSV(dpd,"dpd.csv");

        // sort the distribution
        TreeMap<Double,Double> map = new TreeMap<>();
        for (int i=0; i < probabilities.size(); i++)
            {
                Double key = outcomes.get(i);
                Double value = probabilities.get(i);
                map.put(key,value);
            }

        // setup arrays
        List<Double> peakLocations = new LinkedList<>();
        List<Double> peakAreas = new LinkedList<>();
        
        Set<Map.Entry<Double,Double>> entrySet = map.entrySet();
        List<Map.Entry<Double,Double>> entryList = new LinkedList<>(entrySet);

        Double lastPr = null;
        Double nextPr = null;
        for (int i=0; i < entryList.size(); i++)
            {
                Map.Entry<Double,Double> entry = entryList.get(i);
                Double thisPr = entry.getValue();
                Double outcome = entry.getKey();

                if ( i < entryList.size() - 1 )
                    nextPr = entryList.get(i+1).getValue();
                else
                    nextPr = null;
                //System.out.println(i + ", " + outcome + " Pr = " +thisPr);
                boolean check = false;
                if ( i == 0 && thisPr > nextPr )
                    check = true;
                else if ( i == entryList.size() - 1 && thisPr > lastPr )
                    check = true;
                else if ( i > 0 && i < entryList.size() && thisPr > lastPr && thisPr > nextPr )
                    check = true;

                if ( check )
                    {
                        List<Double> theseProbabilities = new LinkedList<>();
                        List<Double> theseOutcomes = new LinkedList<>();
                        double sum = 0.0;
                        //System.out.println("peak is at " + i + ", x = " + outcome);
                        //System.out.println("range: " + Math.max(0,i-PEAK_SIZE) + " to " + Math.min(entryList.size()-1,i+PEAK_SIZE) );
                        for (int j=Math.max(0,i-PEAK_SIZE); j < Math.min(entryList.size()-1,i+PEAK_SIZE); j++)
                            {
                                Map.Entry<Double,Double> entry2 = entryList.get(j);
                                Double key = entry2.getKey();
                                Double value = entry2.getValue();
                                theseProbabilities.add(value);
                                theseOutcomes.add(key);
                                sum += value;
                            }
                        //System.out.println("area = " + sum);
                        if ( sum > PEAK_THRESHOLD )
                            {
                                DiscreteProbabilityDistribution<Double> newDPD = new DiscreteProbabilityDistribution<Double>(theseOutcomes, theseProbabilities);
                                double expectedValue = getExpectedValue(newDPD);
                                //System.out.println(">>> " + expectedValue);
                                peakLocations.add(expectedValue);
                                peakAreas.add(sum);
                            }
                    }

                lastPr = thisPr;
            }

        // if there's only one peak, simply return the expected value of the distribution
        List<Double> returnList = new LinkedList<>();
        if ( peakLocations.size() <= 1 )
            returnList.add(getExpectedValue(dpd));

        // if there are multiple entries,
        // sort peaks by largest area and return at most MAX_PEAKS peak locations
        else
            {
                TreeMap<Double,Double> peakMap = new TreeMap<>();
                for (int i=0; i < peakLocations.size(); i++)
                    {
                        Double key = peakAreas.get(i);
                        Double value = peakLocations.get(i);
                        peakMap.put(key,value);
                    }
                //System.out.println(peakMap);

                // iterate backwards
                entrySet = peakMap.entrySet();
                entryList = new ArrayList<>(entrySet);
                ListIterator<Map.Entry<Double,Double>> iterator = entryList.listIterator(entryList.size());
                int count = 0;
                while ( iterator.hasPrevious() && count < MAX_PEAKS )
                    {
                        Map.Entry<Double,Double> entry = iterator.previous();
                        returnList.add(entry.getValue());
                        count++;
                    }
            }

        // return the result
        //System.out.println("final: ");
        //System.out.println(returnList);
        return ImmutableList.copyOf(returnList);
    }

    /**
     * Writes a comma separated value file for the given histogram.
     * @param DPD the distribution to be described
     * @param filename the filename to write the CSV to
     */
    public static void writeCSV(DiscreteProbabilityDistribution<Double> DPD, String filename)
    {
        String CSVstring = "";
        List<Double> probabilities = DPD.inputProbabilities;
        List<Double> outcomes = DPD.outcomes;
        for (int i=0; i < probabilities.size(); i++)
            CSVstring += outcomes.get(i) + "," + probabilities.get(i) + "\n";
        InputFileFormat.writeStringToDisk(CSVstring,filename);
    }

    /**
     * Writes a debug file containing the atoms in this rotamer
     * @param r the residue containing this Rotamer
     * @param filename the name of the file to write to
     */
    public void dump(Residue r, String filename)
    {
        List<Atom> dumpAtoms = new LinkedList<>(atoms);
        dumpAtoms.add(r.chis.get(0).atom1);
        dumpAtoms.add(r.chis.get(0).atom2);
        String dumpString = "#\n\ntitle\n\n0 1\n";
        for (Atom a : dumpAtoms)
            dumpString += a.toString() + "\n";
        dumpString += "\n\n";
        InputFileFormat.writeStringToDisk(dumpString, filename);
    }

    /**
     * Determines whether two atoms belonging to two rotamers might be hydrogen bonded.  A hydrogen bond is defined as
     * X...H-X, where X=H,N,S (defined in {@link #HBOND_ELEMENTS}) and X has a free lone pair.  Additionally,
     * the X...H-X angle must be at least {@link #HBONDED_MINIMUM_ANGLE}#.  Nulls are not checked.  We don't check
     * the H-bond distance here, because we assume this will only be called if the distance is appropriate.
     * a1 and a2 are assumed to be in different sidechains and not directly bonded.  These H-bond checks are optimized
     * for the standard amino acids and may fail if there are unusual ones present.
     * @param r1 the rotamer the first atom belongs to
     * @param a1 the first atom
     * @param r2 the rotamer the second atom belongs to
     * @param a2 the second atom
     * @return true if there might be a hydrogen bond
     */
    public static boolean possibleHydrogenBond(Rotamer r1, Atom a1, Rotamer r2, Atom a2)
    {
        // get elements
        Element e1 = a1.element;
        Element e2 = a2.element;

        // check case:
        // X  ... H --- X
        // a1    a2     polarAtom
        if ( HBOND_ELEMENTS.contains(e1) && e2 == Element.HYDROGEN && r2.nextToPolar(a2) != null )
            {
                Atom polarAtom = r2.nextToPolar(a2);

                // check the hydrogen-bond angle
                double angle = Molecule.getAngle(a1, a2, polarAtom);
                if ( angle < HBONDED_MINIMUM_ANGLE )
                    return false;

                // check a1 is the correct valence
                if ( e1 == Element.NITROGEN && r1.numberOfNeighbors(a1) > 2 )
                    return false;

                // this must be a hydrogen bond
                return true;
            }

        // check case:
        // X  ... H  --- X
        // a2     a1     polarAtom
        else if ( HBOND_ELEMENTS.contains(e2) && e1 == Element.HYDROGEN && r1.nextToPolar(a1) != null )
            {
                Atom polarAtom = r1.nextToPolar(a1);

                // check the hydrogen-bond angle
                double angle = Molecule.getAngle(a2,a1,polarAtom);
                if ( angle < HBONDED_MINIMUM_ANGLE )
                    return false;

                // check a2 is the correct valence
                if ( e2 == Element.NITROGEN && r2.numberOfNeighbors(a2) > 2 )
                    return false;

                // this must be a hydrogen bond
                return true;
            }
        return false;
    }

    /**
     * Tells us if atom a is next to a polar atom.  For determining whether there is a hydrogen bond.
     * @param a an atom in this sidechain
     * @return null if this isn't next to a polar atom, otherwise the polar atom
     */
    public Atom nextToPolar(Atom a)
    {
        int index = atoms.indexOf(a);
        return nextToPolar.get(index);
    }

    /**
     * Tells us how many bonded neighbors atom a has.  For determining whether there is a hydrogen bond.
     * @param a an atom in this sidechain
     * @return the number of bonded neighbors
     */
    public int numberOfNeighbors(Atom a)
    {
        int index = atoms.indexOf(a);
        return numberOfNeighbors.get(index);
    }

    /**
     * Returns the interesting atom in this Rotamer.
     * Only the first hit is returned.  Intended to be used on structures that are not using the close
     * contact atom types.  See atom typing.cdx.
     * @return the interesting atom if present, null otherwise
     */
    public Atom getInterestingAtom()
    {
        for (Atom a : atoms)
            {
                if ( a.tinkerAtomType == 126 || a.tinkerAtomType == 130 ||   // histidine pi nitrogens
                     a.tinkerAtomType == 401                               ) // methanol HO
                    return a;
            }
        return null;
    }

    @Override
    public boolean equals(Object obj)
    {
        if ( obj == null )
            return false;
        if ( obj == this )
            return true;
        if ( !(obj instanceof Rotamer) )
            return false;

        Rotamer r = (Rotamer)obj;
        if ( sequenceIndex == r.sequenceIndex && 
             Objects.equals(atoms, r.atoms)      )
             /*Objects.equals(chis, r.chis) &&
             Objects.equals(protoAminoAcid, r.protoAminoAcid) &&
             Objects.equals(nextToPolar, r.nextToPolar) &&
             Objects.equals(numberOfNeighbors, r.numberOfNeighbors) &&*/
            return true;
        return false;
    }

    @Override
    public int hashCode()
    {
        //return Objects.hash(atoms, sequenceIndex, chis, protoAminoAcid, nextToPolar, numberOfNeighbors);
        return Objects.hash(sequenceIndex, atoms);
    }

    @Override
    public String toString()
    {
        return String.format("[%d] %s: %s", sequenceIndex, protoAminoAcid.molecule.name, chis.toString());
        //for (Atom a : atoms)
        //    returnString += a.toString() + "\n";
    }
}
