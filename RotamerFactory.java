import java.util.*;
import java.io.*;
import org.apache.commons.math3.geometry.euclidean.threed.*;
import org.jgrapht.*;
import org.jgrapht.graph.*;
import org.jgrapht.alg.*;
import com.google.common.collect.*;
import java.util.concurrent.*;

/**
 * This class makes rotamers from a given peptide.
 */
public class RotamerFactory
{
    /** This class is not instantiable. */
    private RotamerFactory()
    {
        throw new IllegalArgumentException("not instantiable");
    }

    /**
     * Creates a single rotamer using the atoms that are already in the residue.
     * @return a list that only contains the current rotamer
     */
    public static List<Rotamer> getOneRotamer(Peptide peptide, Residue residue, boolean includeHN)
    {
        // return one rotamer that is just the current rotamer
        List<Atom> singleRotamerAtoms = new ArrayList<>(getSidechainAtoms(peptide, residue, includeHN));
        
        // determine the chis
        List<Double> singleRotamerChiList = null;
        if ( residue.chis.size() == 0 )
            singleRotamerChiList = ImmutableList.of();
        else
            {
                singleRotamerChiList = new ArrayList<>(residue.chis.size());
                for (ProtoTorsion t : residue.chis)
                    singleRotamerChiList.add(t.getDihedralAngle());
                singleRotamerChiList = ImmutableList.copyOf(singleRotamerChiList);
            }

        int sequenceIndex = peptide.sequence.indexOf(residue);
        if ( sequenceIndex == -1 )
            throw new IllegalArgumentException("residue not found in sequence");
        Rotamer singleRotamer = new Rotamer(singleRotamerAtoms, sequenceIndex, singleRotamerChiList, residue.description);
        List<Rotamer> returnList = new ArrayList<>(1);
        returnList.add(singleRotamer);
        return returnList;
    }

    /**
     * Finds the possible rotamers.  The method will draw chis for normal residues (not TS, histidine or hairpin)
     * Assumes that the residue is already of the correct type.  This method does not check the resulting rotamers for
     * clashes with anything (backbone or other rotamers).
     * If there are no rotamers possible, an empty list is returned.
     * @param peptide the peptide to analyze
     * @param residue the residue to generate rotamers for
     * @param includeHN if set to true, the backboneHN of this residue will be included in the atoms list of the resulting rotamers
     * @param inputChis if non-null, these are lists of chi values in degrees that will be used
     * @return a list of rotamers that are possible at the specified position
     */
    public static List<Rotamer> generateRotamers(Peptide peptide, Residue residue, boolean includeHN, List<List<Double>> inputChis)
    {
        // check this is an appropriate residue
        if ( !peptide.sequence.contains(residue) )
            throw new IllegalArgumentException("residue not in peptide");

        // setup some variables
        int sequenceIndex = peptide.sequence.indexOf(residue);
        String description = residue.description;
        AminoAcid aminoAcid = residue.aminoAcid;
        if ( inputChis == null && aminoAcid.chirality == Chirality.D )
            throw new IllegalArgumentException("can't use this method to generate rotamers for " + description);
        List<Rotamer> returnList = new ArrayList<>(); // result rotamers will be placed here

        // get possible chi angles
        List<List<Double>> possibleChis = null;
        if ( inputChis != null )
            {
                // no chis so return no rotamers
                if ( inputChis.size() == 0 )
                    return returnList;

                // the method has been called with specific chi values, so try to use those
                // first check if there's the expected number of chis
                for (List<Double> list : inputChis)
                    if ( list.size() != residue.chis.size() )
                        throw new IllegalArgumentException("unexpected number of chis");
                possibleChis = inputChis;
            }
        else
            {
                // the method has not been called with specific chi values
                if ( residue.aminoAcid.rotamerType == AminoAcid.RotamerType.HAS_NO_ROTAMERS )
                    {
                        // return one rotamer that is just the current rotamer
                        return getOneRotamer(peptide, residue, includeHN);
                    }
                // D amino acids are not supported 
                else if ( residue.aminoAcid.chirality == Chirality.D )
                    throw new IllegalArgumentException("no data for D amino acids");
                // draw normal chis
                possibleChis = RotamerSummarizer.getPossibleRotamers(residue);
            }
        if ( possibleChis.size() == 0 )
            throw new IllegalArgumentException("no chis for " + description);

        // get the torsions that describe where the sidechain is
        LinkedList<ProtoTorsion> oldProtoTorsions = new LinkedList<>(residue.chis);

        // get the connectivity graph
        // if this is a proline, make a copy of the graph and disconnect the chi3 bond to the backbone
        SimpleWeightedGraph<Atom,DefaultWeightedEdge> connectivity = peptide.connectivity;
        if ( aminoAcid.isProline() )
            {
                // clone the graph
                connectivity = Molecule.cloneConnectivity(connectivity);

                // remove the bond that connects the proline ring to the backbone at the beta carbon
                ProtoTorsion chi3 = residue.chis.get(2);
                Atom atom3 = chi3.atom3;
                Atom atom4 = chi3.atom4;
                DefaultWeightedEdge e = connectivity.removeEdge(atom3,atom4);
                if ( e == null )
                    throw new NullPointerException("unexpected null edge for proline connection");
            
                // don't set the last chi angle
                oldProtoTorsions.removeLast();
            }
        
        // get the atoms in the sidechain
        Pair<Atom,Atom> prochiralConnection = residue.prochiralConnection;
        List<Atom> allAtoms = new ArrayList<>(getHalfGraph(connectivity, prochiralConnection.getFirst(), prochiralConnection.getSecond()));
        if ( includeHN && residue.HN != null )
            allAtoms.add(residue.HN);

        // add the backbone atoms we need for the first torsions
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
        for (int j=0; j < possibleChis.size(); j++)
            {
                List<Double> requestedAngles = possibleChis.get(j);

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

                // return the new Rotamer
                Rotamer newRotamer = new Rotamer(newAtoms, sequenceIndex, requestedAngles, description);
                returnList.add(newRotamer);
            }

        // return the result
        return returnList;
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
     * Given a bond between includeAtom and excludeAtom, returns a set of atoms
     * containing all the atoms on the includeAtom side of the bond, including
     * includeAtom.  Will throw an exception if includeAtom and excludeAtom are not
     * directly bonded.  Will throw an exception if includeAtom and excludeAtom
     * form a ring.
     * @param connectivity the connectivity graph
     * @param excludeAtom this atom will not be included in the result
     * @param includeAtom this atom will be included in the result
     * @return the atoms on the includeAtom side of the graph
     */
    public static Set<Atom> getHalfGraph(SimpleWeightedGraph<Atom,DefaultWeightedEdge> connectivity, Atom excludeAtom, Atom includeAtom)
    {
        // if these atoms are not directly bonded, then return an empty set
        if ( !connectivity.containsEdge(excludeAtom,includeAtom) )
            throw new IllegalArgumentException("atoms not connected"); 

        // preform a breadth-first search of one branch of the graph only
        Set<Atom> returnSet = new HashSet<Atom>();
        LinkedList<Atom> searchQueue = new LinkedList<Atom>();
        Set<Atom> searched = new HashSet<>();
        for (Atom a : getAdjacentAtoms(connectivity,includeAtom))
            {
                searchQueue.add(a);
                returnSet.add(a);
            }
        searchQueue.remove(excludeAtom);
        returnSet.remove(excludeAtom);
        returnSet.add(includeAtom);
        searched.add(includeAtom);
        
        Atom lastNode = includeAtom;
        while (searchQueue.size() > 0)
            {
                Atom currentNode = searchQueue.remove();
                Set<Atom> adjacent = getAdjacentAtoms(connectivity,currentNode);
                adjacent.remove(lastNode);
                for (Atom a : adjacent)
                    {
                        if ( a == excludeAtom)
                            {
                                // if the excluded atom is found, this is a ring!
                                throw new IllegalArgumentException("RotamerFactory.getHalfGraph found a ring");
                            }

                        // if this is an atom we haven't already searched, mark it
                        if ( ! searched.contains(a) )
                            {
                                searchQueue.add(a);
                                searched.add(a);
                                returnSet.add(a);
                            }
                    }
            }
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

    /**
     * Returns the sidechain atoms for a given residue.
     * @param peptide the peptide the residue is in
     * @param residue the residue to get the sidechain atoms of
     * @param includeHN if true, the backbone HN will be included in the result
     * @return the sidechain atoms of the specified residue
     */
    public static Set<Atom> getSidechainAtoms(Peptide peptide, Residue residue, boolean includeHN)
    {
        Set<Atom> sidechainAtoms = null;
        Atom atomA = residue.prochiralConnection.getFirst();
        Atom atomB = residue.prochiralConnection.getSecond();
        if ( residue.aminoAcid.isProline() )
            {
                // this is a proline, so create another molecule where the proline
                // sidechain has been disconnected
                ProtoTorsion chi3 = residue.chis.get(2);
                Atom atom3 = chi3.atom3;
                Atom atom4 = chi3.atom4;

                // temporarily disconnect proline ring bond
                SimpleWeightedGraph<Atom,DefaultWeightedEdge> connectivity = Molecule.cloneConnectivity(peptide.connectivity);
                DefaultWeightedEdge edge = connectivity.removeEdge(atom3, atom4);
                if ( edge == null )
                    throw new NullPointerException("expected edge for proline");

                // traverse the sidechain subgraph
                sidechainAtoms = getHalfGraph(connectivity, atomA, atomB);
            }
        else
            {
                // just add the sidechain atoms
                sidechainAtoms = peptide.getHalfGraph(atomA, atomB);
                if ( includeHN )
                    {
                        if ( residue.HN == null )
                            throw new NullPointerException("expected HN");
                        sidechainAtoms.add(residue.HN);
                    }
            }
        return sidechainAtoms;
    }

    /**
     * Returns the atoms that are in the backbone.  Backbone atoms are defined as those that will never move during
     * rotamer packing.  If a position has no rotamers, all atoms at the position are considered part of the backbone.
     * If a position is marked as variable, its sidechain atoms will not be included in the result.
     * @param peptide the peptide to analyze
     * @param variablePositions indices in the peptide sequence where we want to vary the rotamers
     * @param includeHN if true, the backbone HNs are considered part of the sidechain
     * @return the atoms in the backbone
     */
    public static List<Atom> getBackboneAtoms(Peptide peptide, List<Integer> variablePositions, boolean includeHN)
    {
        List<Atom> backboneAtoms = new ArrayList<>(peptide.contents);
        for (int i=0; i < peptide.sequence.size(); i++)
            {
                Residue residue = peptide.sequence.get(i);
                if ( variablePositions.contains(i) && residue.aminoAcid.rotamerType != AminoAcid.RotamerType.HAS_NO_ROTAMERS )
                    backboneAtoms.removeAll( getSidechainAtoms(peptide, residue, includeHN) );
            }
        return backboneAtoms;
    }

    /** how many points to use to search for TS rotamers */
    public static final int TS_GRID_SIZE = 13;

    /**
     * Finds the chi values for a transition state amino acid, under the constraint that the
     * anionic oxygen must be near a backbone HN.  Clashes with the backbone are crudely checked.
     * Assumes this is a beta hairpin.
     * @param peptide the peptide to analyze
     * @param residue the residue to mutate
     * @param includeHN if true, the backbone HN will be included in the rotamer
     * @return lists of lists (chi1, chi2) for transition rotamers that satisfy the constraint
     */
    public static List<List<Double>> getTransitionStateChis(Peptide peptide, Residue residue, boolean includeHN)
    {
        // check invariants
        if ( residue.description.indexOf("transition_state") == -1 )
            throw new IllegalArgumentException("must be a transition state");
        int sequenceIndex = peptide.sequence.indexOf(residue);
        if ( sequenceIndex == -1 )
            throw new IllegalArgumentException("bad sequence index");
 
        // get the backbone atoms
        int sequenceLength = peptide.sequence.size();
        if ( sequenceLength % 2 != 0 )
            throw new IllegalArgumentException("expecting even sequence length");
        List<Integer> variablePositions = new ArrayList<>(sequenceLength-2);
        int forbiddenIndex = (sequenceLength/2) - 1;
        if ( sequenceIndex == forbiddenIndex || sequenceIndex == forbiddenIndex + 1 )
            throw new IllegalArgumentException("can't place TS at a hairpin position");
        for (int i=0; i < sequenceLength; i++)
            {
                if ( i == forbiddenIndex )
                    {
                        i++;
                        continue;
                    }
                variablePositions.add(i);
            }
        List<Atom> backboneAtoms = getBackboneAtoms(peptide, variablePositions, includeHN);
        //for (Atom a : backboneAtoms)
        //    System.out.print(peptide.getAtomString(a) + ", " );
        //System.out.println();

        // get the rotation axes
        ProtoTorsion chi1torsion = residue.chis.get(0);
        ProtoTorsion chi2torsion = residue.chis.get(1);

        // get the sidechain positions
        Pair<Atom,Atom> prochiralConnection = residue.prochiralConnection;
        List<Atom> allAtoms = new ArrayList<>(peptide.getHalfGraph(prochiralConnection.getFirst(), prochiralConnection.getSecond()));
        allAtoms.add(chi1torsion.atom1);
        allAtoms.add(chi1torsion.atom2);

        // add the backboneHN if requested
        if ( residue.HN == null )
            throw new NullPointerException("null HN for transition state");
        if ( includeHN )
            allAtoms.add(residue.HN);

        // we don't want to check atoms that are part of the backbone for clashes
        // but just in case there are numerical issues we refer to them by index instead of reference
        HashSet<Integer> ignoreAtomIndices = new HashSet<>();
        ignoreAtomIndices.add(allAtoms.indexOf(chi1torsion.atom1));
        ignoreAtomIndices.add(allAtoms.indexOf(chi1torsion.atom2));
        ignoreAtomIndices.add(allAtoms.indexOf(chi1torsion.atom3));
        if ( includeHN && residue.HN != null )
            ignoreAtomIndices.add(allAtoms.indexOf(residue.HN));

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
                if ( a.type1 == 408 )
                    {
                        if ( atomOindex != null )
                            throw new IllegalArgumentException("two TS oxygens found");
                        atomOindex = allAtoms.indexOf(a);
                    }
            }
        if ( atomOindex == null )
            throw new NullPointerException("null oxygen index");

        // find compatible HNs
        // this must be the HN adjacent to the position and the HN must not be hydrogen
        // bonded in the sheet; this means that there may be no solution for some residues
        List<Integer> HNindices = new ArrayList<>(); // sequence indices of adjacent residues with accessible, non-sheet HNs
        for (int k=sequenceIndex-1; k < sequenceIndex+2; k++)
            {
                if ( k >= 0 && k < sequenceLength &&                                     // must be inside the sequence
                     k != forbiddenIndex && k != forbiddenIndex+1 )                      // must not be a hairpin position
                    {
                        if ( ( sequenceLength / 2 ) % 2 == 0 )                           // n/2 even; n = 8, 12, 16, ...
                            {
                                if ( k < forbiddenIndex && k % 2 == 1 )                  // before hairpin, odd
                                    HNindices.add(k);
                                else if ( k > forbiddenIndex + 1 && k % 2 == 0 )         // after hairpin, even
                                    HNindices.add(k);
                            }
                        else if ( ( sequenceLength / 2 ) % 2 == 1 )                      // n/2 odd; n = 10, 14, 18, ...
                            {
                                if ( k < forbiddenIndex && k % 2 == 0 )                  // before hairpin, even
                                    HNindices.add(k);
                                else if ( k > forbiddenIndex + 1 && k % 2 == 1 )         // after hairpin, odd
                                    HNindices.add(k);
                            }
                    }
            }
        //System.out.println(sequenceIndex + " : " + HNindices);

        // rotate chi1 and chi2 on a grid
        // check for the desired contact
        List<List<Double>> returnList = new ArrayList<>();
        double stepSize = 360.0 / (TS_GRID_SIZE-1.0);
        double chi1 = -180.0;
        for (int i=0; i < TS_GRID_SIZE; i++)
            {
                if ( HNindices.size() == 0 )
                    break;
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
                        // check for N-H...O contact distance and angle
                        boolean acceptable = false;
                        for (Integer HNindex : HNindices)
                            {
                                Residue peptideResidue = peptide.sequence.get(HNindex);
                                Atom atomN = peptideResidue.N;
                                Atom atomHN = peptideResidue.HN;

                                // check distance and angle
                                double distance = Vector3D.distance(atomO.position, atomHN.position);
                                if ( distance < Settings.MAXIMUM_HBOND_DISTANCE && distance > 1.25 &&
                                     Molecule.getAngle(atomO.position, atomHN.position, atomN.position) > Settings.MINIMUM_HBOND_ANGLE )
                                    {
                                        acceptable = true;
                                        break;
                                    }
                            }

                        // if acceptable, check for clashes and add it to the list of chis to return if it's ok
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
                                                if ( distance < Settings.MINIMUM_INTERATOMIC_DISTANCE )
                                                    {
                                                        //System.out.println(rotamerAtom.type1);
                                                        //System.out.println(peptide.getAtomString(backboneAtom) + " " + distance);
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
     * Determines if the specified residue is up or down.  Exceptions will be thrown if this isn't a hairpin or a hairpin position is given.
     * @param sequenceLength the length of the sequence
     * @param sequenceIndex the index of the residue
     * @return true if the substituent is up, false if it is down (the absolute direction is arbitrary but is always consistent)
     */
    public static boolean isUp(int sequenceLength, int sequenceIndex)
    {
        if ( sequenceLength % 2 != 0 )
            throw new IllegalArgumentException("expecting even sequence length");
        int halfway = sequenceLength/2 - 1; // the index of the first hairpin position (i.e., d-pro)
        if ( sequenceIndex == halfway || sequenceIndex == halfway + 1 )
            throw new IllegalArgumentException("this is a hairpin position");
        if ( sequenceIndex < 0 || sequenceIndex > sequenceLength-1 )
            throw new IllegalArgumentException("sequence index exceeds sequence length");
        if ( (sequenceLength/2) % 2 == 0 )
            {
                // n = 12, 16, ...
                if ( sequenceIndex < halfway )
                    {
                        if ( sequenceIndex % 2 == 0 )
                            return false;
                        return true;
                    }
                else
                    {
                        if ( sequenceIndex % 2 == 0 )
                            return true;
                        return false;
                    }
            }
        else if ( (sequenceLength/2) % 2 == 1 )
            {
                // n = 10, 14, ...
                if ( sequenceIndex < halfway )
                    {
                        if ( sequenceIndex % 2 == 0 )
                            return true;
                        return false;
                    }
                else
                    {
                        if ( sequenceIndex % 2 == 0 )
                            return false;
                        return true;
                    }
 
            }
        throw new IllegalArgumentException("unreachable");
    }

    /**
     * Takes a poly-gly template and creates interesting pairs of transition states and histidines.
     * @param betaSheet the poly-gly template
     * @param includeHN whether to include the backboneHN in the sidechain
     * @return pairs of TS and histidine rotamers
     */
    public static List<Pair<Rotamer,Rotamer>> generateInterestingPairs(Peptide betaSheet, boolean includeHN)
    {
        // create peptides that have each kind of transition state at each non-hairpin position
        List<ProtoAminoAcid> tsTemplates = new ArrayList<>(3);
        tsTemplates.add(ProtoAminoAcidDatabase.getTemplate("ts1"));
        tsTemplates.add(ProtoAminoAcidDatabase.getTemplate("ts2"));
        tsTemplates.add(ProtoAminoAcidDatabase.getTemplate("ts3"));
        
        int sequenceLength = betaSheet.sequence.size();
        if ( sequenceLength % 2 != 0 )
            throw new IllegalArgumentException("even sequence length expected");
        int forbiddenIndex = (sequenceLength/2) - 1;
        Map<Peptide,Integer> transitionStatePeptides = new HashMap<>(); // peptides mapped to TS sequence index
        for (int i=0; i < sequenceLength; i++)
            {
                if ( i == forbiddenIndex )
                    {
                        i++;
                        continue;
                    }
                for (ProtoAminoAcid tsTemplate : tsTemplates)
                    {
                        Peptide newPeptide = SidechainMutator.mutateSidechain(betaSheet, betaSheet.sequence.get(i), tsTemplate);
                        transitionStatePeptides.put(newPeptide, i);
                    }
            }

        // generate pairs of transition state rotamers and histidine rotamers
        List<Pair<Rotamer,Rotamer>> returnList = new ArrayList<>();
        for (Peptide tsPeptide : transitionStatePeptides.keySet())
            {
                // generate transition state rotamers
                int TSindex = transitionStatePeptides.get(tsPeptide);
                List<List<Double>> TSchis = getTransitionStateChis(tsPeptide, tsPeptide.sequence.get(TSindex), includeHN);
                List<Rotamer> transitionStateRotamers = generateRotamers(tsPeptide, tsPeptide.sequence.get(TSindex), includeHN, TSchis);
                //System.out.printf("%d TS rotamers generated\n", transitionStateRotamers.size());

                // generate histidine rotamers
                for (Rotamer TSrotamer : transitionStateRotamers)
                    {
                        List<Pair<Rotamer,Rotamer>> rotamerPairs = getHistidineRotamerPairs(tsPeptide, TSrotamer, includeHN);
                        //System.out.printf("   %d histidine rotamers generated\n", rotamerPairs.size());
                        returnList.addAll(rotamerPairs);
                    }
            }
        return returnList;
    }

    /** how many points to use to search for histidine rotamers */
    public static final int HISTIDINE_GRID_SIZE = 13;

    /**
     * Finds the chi values for a histidine amino acid, under the constraint that the pi-nitrogen must be near a
     * transition state hydroxyl.
     * @param inputPeptide the peptide to generate histidine rotamers for
     * @param transitionStateRotamer the transition state rotamer to find compatible histidines with
     * @param includeHN if true, include the HN in the resulting rotamers
     * @return pairs of transition state and histidine rotamers that should be reactive
     */
    public static List<Pair<Rotamer,Rotamer>> getHistidineRotamerPairs(Peptide inputPeptide, Rotamer transitionStateRotamer, boolean includeHN)
    {
        // make a list of where we can put histidines
        // the histidine should be on the same face of the hairpin as the transition state
        // and can't be where the transition state or hairpin residues are
        int sequenceLength = inputPeptide.sequence.size();
        int sequenceIndex = transitionStateRotamer.sequenceIndex;
        boolean transitionStateIsUp = isUp(sequenceLength, sequenceIndex); // figure out whether the transition state rotamer is up or down
        int forbiddenIndex = (sequenceLength/2) - 1;
        List<Integer> availablePositions = new ArrayList<>(sequenceLength/2);
        for (int i=0; i < sequenceLength; i++)
            {
                // skip over hairpin positions
                if ( i == forbiddenIndex )
                    {
                        i++;
                        continue;
                    }
                // can't place a histidine at the same position as the transition state
                if ( i == sequenceIndex )
                    continue;

                // if this residue is on the same side as the transition state, mark the position as open
                boolean positionIsUp = isUp(sequenceLength,i);
                if ( ( transitionStateIsUp && positionIsUp ) || ( ! transitionStateIsUp && ! positionIsUp ) )
                    availablePositions.add(i);
            }
        if ( availablePositions.size() == 0 )
            throw new IllegalArgumentException("should have at least one available position");

        // make peptides containing a histidine at each valid position, noting that we need both HID and HIE tautomers
        // map is from peptides to the sequence index of the histidine
        ProtoAminoAcid histidine_HID = ProtoAminoAcidDatabase.getTemplate("histidine_hd");
        ProtoAminoAcid histidine_HIE = ProtoAminoAcidDatabase.getTemplate("histidine_he");
        Map<Peptide,Integer> histidinePeptides = new HashMap<>();
        for (Integer i : availablePositions)
            {
                histidinePeptides.put(SidechainMutator.mutateSidechain(inputPeptide, inputPeptide.sequence.get(i), histidine_HID), i);
                histidinePeptides.put(SidechainMutator.mutateSidechain(inputPeptide, inputPeptide.sequence.get(i), histidine_HIE), i);
            }

        // get the transition state hydroxyl atoms
        Atom atomH = locateSingleAtom(ImmutableSet.of(401), transitionStateRotamer.atoms);
        Atom atomO = locateSingleAtom(ImmutableSet.of(402), transitionStateRotamer.atoms);

        // for each peptide, rotate the sidechain on a grid and see if it gets close to the transition state hydroxyl
        List<Pair<Rotamer,Rotamer>> returnList = new ArrayList<>();
        for (Peptide histidinePeptide : histidinePeptides.keySet())
            {
                // get backbone atoms
                // the backbone will include everything but the histidine sidechain
                // this makes sense if the peptide is a poly-gly with one transition state
                // if we run into glycine Hs, we will have a problem with any other sidechain anyways
                int histidineResidueIndex = histidinePeptides.get(histidinePeptide);
                List<Atom> backboneAtoms = getBackboneAtoms(histidinePeptide, ImmutableList.of(histidineResidueIndex), includeHN);

                // get the rotation axes
                Residue residue = histidinePeptide.sequence.get(histidineResidueIndex);
                ProtoTorsion chi1torsion = residue.chis.get(0);
                ProtoTorsion chi2torsion = residue.chis.get(1);

                // get the sidechain positions
                Pair<Atom,Atom> prochiralConnection = residue.prochiralConnection;
                List<Atom> allAtoms = new ArrayList<>(histidinePeptide.getHalfGraph(prochiralConnection.getFirst(), prochiralConnection.getSecond()));
                allAtoms.add(chi1torsion.atom1);
                allAtoms.add(chi1torsion.atom2);

                // add the backboneHN if requested
                if ( residue.HN == null )
                    throw new NullPointerException("null HN for histidine");
                if ( includeHN )
                    allAtoms.add(residue.HN);

                // we don't want to check atoms that are part of the backbone for clashes
                // but just in case there are numerical issues we refer to them by index instead of reference
                HashSet<Integer> ignoreAtomIndices = new HashSet<>();
                ignoreAtomIndices.add(allAtoms.indexOf(chi1torsion.atom1));
                ignoreAtomIndices.add(allAtoms.indexOf(chi1torsion.atom2));
                ignoreAtomIndices.add(allAtoms.indexOf(chi1torsion.atom3));
                if ( includeHN && residue.HN != null )
                    ignoreAtomIndices.add(allAtoms.indexOf(residue.HN));

                // the atoms to rotate for each torsion
                List<Atom> chi1atoms = new ArrayList<>(histidinePeptide.getHalfGraph(chi1torsion.atom2, chi1torsion.atom3));
                List<Atom> chi2atoms = new ArrayList<>(histidinePeptide.getHalfGraph(chi2torsion.atom2, chi2torsion.atom3));
            
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

                // get the index of the histidine pi nitrogen
                Atom histidinePiNitrogen = locateSingleAtom(ImmutableSet.of(130,126), allAtoms); 
                Integer histidinePiNitrogenIndex = allAtoms.indexOf(histidinePiNitrogen);

                // rotate chi1 and chi2 on a grid and check for the desired contact
                double stepSize = 360.0 / (HISTIDINE_GRID_SIZE-1.0);
                double chi1 = -180.0;
                for (int i=0; i < HISTIDINE_GRID_SIZE; i++)
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
                        for (int j=0; j < HISTIDINE_GRID_SIZE; j++)
                            {
                                // try this rotamer
                                List<Atom> thisAtoms2 = setDihedral(thisAtoms, chi2atoms, chi2torsion, chi2);

                                // get the atoms in the putative hydrogen bond
                                Atom atomN = thisAtoms2.get(histidinePiNitrogenIndex);

                                // check that the oxygen forms a hydrogen bond to a backbone NH
                                // check for N-H...O contact distance and angle
                                double distance = Vector3D.distance(atomN.position, atomH.position);
                                if ( distance < Settings.MAXIMUM_HBOND_DISTANCE && distance > 1.25 &&
                                     Molecule.getAngle(atomN.position, atomH.position, atomO.position) > 150.0 )
                                     //Molecule.getAngle(atomN.position, atomH.position, atomO.position) > Settings.MINIMUM_HBOND_ANGLE )
                                    {
                                        // check for backbone atom clashes
                                        HashSet<Atom> ignoreAtoms = new HashSet<>();
                                        for (Integer index : ignoreAtomIndices)
                                            ignoreAtoms.add(thisAtoms2.get(index));

                                        // check if this rotamer now clashes with the backbone
                                        boolean acceptable = true;
                                        checking:
                                        for (Atom rotamerAtom : thisAtoms2)
                                            {
                                                if ( ignoreAtoms.contains(rotamerAtom) )
                                                    continue;

                                                for (Atom backboneAtom : backboneAtoms)
                                                    {
                                                        if ( Vector3D.distance(rotamerAtom.position, backboneAtom.position) < Settings.MINIMUM_INTERATOMIC_DISTANCE )
                                                            {
                                                                acceptable = false;
                                                                break checking;
                                                            }
                                                    }
                                            }
                                        if ( acceptable )
                                            {
                                                // build the actual rotamer
                                                List<List<Double>> thisChis = new ArrayList<>(1);
                                                thisChis.add(ImmutableList.of(chi1, chi2));
                                                Rotamer histidineRotamer = generateRotamers(histidinePeptide, residue, includeHN, thisChis).get(0);
                                                Pair<Rotamer,Rotamer> thisRotamerPair = new Pair<>(transitionStateRotamer, histidineRotamer);
                                                returnList.add(thisRotamerPair);
                                            }
                                    }

                                // increment chi2
                                chi2 += stepSize;
                            }
                        
                        // increment chi1
                        chi1 += stepSize;
                    }
            }

        // return pairs of interesting rotamers: (transition state rotamer, histidine rotamer)
        return returnList;
    }

    /**
     * Trys to add arginine to the specified sheets.
     * For every peptide specified, it tries to place an arginine at every non-hairpin, non-transition-state location.
     * If the rotamer doesn't clash and contains a hydrogen bond with the transition state it is minimized.  All operations
     * in this method occur in serial.
     * @param sheets the input peptides, which should be on the close contact forcefield, and have one TS and one histidine
     * @return peptides with arg on the close contact forcefield
     */
    public static List<Peptide> addArginine(List<Peptide> sheets)
    {
        List<Peptide> returnList = new ArrayList<>();
        ProtoAminoAcid arginine = ProtoAminoAcidDatabase.getTemplate("arginine");
        for (Peptide sheetPeptide : sheets)
            {
                // determine what positions are valid for placing arginine
                // will only place this on the same side as the transition state
                Peptide inputPeptide = HydrogenBondMutator.unmutate(sheetPeptide);
                List<Integer> validPositions = new ArrayList<>();
                int sequenceLength = inputPeptide.sequence.size();
                int forbiddenIndex = (sequenceLength/2) - 1;
                Integer TSindex = null;
                Atom transitionStateOatom = null;
                for (int i=0; i < sequenceLength; i++)
                    {
                        Residue residue = inputPeptide.sequence.get(i);
                        if ( residue.aminoAcid == AminoAcid.TS )
                            {
                                if ( TSindex == null )
                                    {
                                        TSindex = i;
                                        for (Atom a : residue.atoms)
                                            {
                                                if ( a.type1 == 408 )
                                                    {
                                                        if ( transitionStateOatom == null )
                                                            transitionStateOatom = a;
                                                        else
                                                            throw new IllegalArgumentException("duplicate TS O atoms");
                                                    }
                                            }
                                    }
                                else
                                    throw new IllegalArgumentException("two TSs not allowed");
                            }
                    }
                if ( TSindex == null )
                    throw new NullPointerException("TS index not found");
                if ( transitionStateOatom == null )
                    throw new NullPointerException("TS O atom not found");
                
                // whether the arginine should be up or down
                // ensures that the arginine will be on the same face as the TS/his pair
                boolean shouldBeUp = isUp(sequenceLength, TSindex);

                for (int i=0; i < sequenceLength; i++)
                    {
                        Residue residue = inputPeptide.sequence.get(i);
                        if ( i == forbiddenIndex )
                            {
                                i++;
                                continue;
                            }
                        if ( residue.aminoAcid == AminoAcid.GLY && isUp(sequenceLength,i) == shouldBeUp )
                            validPositions.add(i);
                    }

                // for each position, introduce an arginine residue
                //System.out.println(validPositions);
                for (Integer index : validPositions)
                    {
                        // mutate to arginine
                        Peptide candidatePeptide = SidechainMutator.mutateSidechain(inputPeptide, inputPeptide.sequence.get(index), arginine);

                        // draw rotamers for this position
                        Residue residue = candidatePeptide.sequence.get(index);
                        List<Rotamer> rotamers = generateRotamers(candidatePeptide, residue, true, null);
                        //System.out.println(rotamers.size() + " rotamers generated");

                        // get the atoms at all the other positions
                        List<Atom> otherAtoms = new ArrayList<>();
                        for (int i=0; i < sequenceLength; i++)
                            {
                                if ( i == index )
                                    continue;
                                otherAtoms.addAll(candidatePeptide.sequence.get(i).atoms);
                            }

                        // for each rotamer, check if it is interesting
                        List<Rotamer> interestingRotamers = new ArrayList<>();
                        for (Rotamer rotamer : rotamers)
                            {
                                boolean interesting = false;
                                for (Atom a : rotamer.atoms)
                                    {
                                        if ( ( a.type1 == 209 || a.type1 == 212 ) && Molecule.getDistance(a, transitionStateOatom) < 2.5 )
                                            {
                                                interesting = true;
                                                break;
                                            }
                                    }
                                if ( !interesting )
                                    {
                                        //System.out.println("not interesting");
                                        continue;
                                    }

                                // if it is interesting, check to see it clashes with any other residue
                                boolean clashes = false;
                                clashCheck:
                                for (Atom a1 : rotamer.atoms)
                                    {
                                        for (Atom a2 : otherAtoms)
                                            {
                                                //if ( Molecule.getDistance(a1,a2) < Settings.MINIMUM_INTERATOMIC_DISTANCE )
                                                if ( Molecule.getDistance(a1,a2) < 1.50 )
                                                    {
                                                        clashes = true;
                                                        //System.out.printf("clash found between %d and %s\n", a1.type1, a2.type1);
                                                        break clashCheck;
                                                    }
                                            }
                                    }

                                // add the rotamer if it's interesting
                                if ( !clashes )
                                    interestingRotamers.add(rotamer);
                            }

                        // add to the final pile of results
                        for (Rotamer rotamer : interestingRotamers)
                            {
                                Peptide reconstituted = Rotamer.reconstitute(candidatePeptide, rotamer);
                                reconstituted = HydrogenBondMutator.mutate(reconstituted);
                                returnList.add(reconstituted);
                            }
                    }
            }
        return returnList;
    }

    /**
     * Places transistion state serine, histidine, and arginine in a beta sheet.
     */
    public static class InterestingJob implements WorkUnit
    {
        public final Peptide peptide;
        public final List<Peptide> interestingPeptides;

        public InterestingJob(Peptide peptide, List<Peptide> interestingPeptides)
        {
            this.peptide = peptide;
            this.interestingPeptides = interestingPeptides;
        }

        public Result call()
        {
            // form interesting tuples
            List<Peptide> interestingPairPeptides = new ArrayList<>();
            List<Pair<Rotamer,Rotamer>> rotamerPairs = generateInterestingPairs(peptide, true);
            for (Pair<Rotamer,Rotamer> rotamerPair : rotamerPairs)
                {
                    Rotamer TSrotamer = rotamerPair.getFirst();
                    Rotamer histidineRotamer = rotamerPair.getSecond();
                    Peptide newPeptide = Rotamer.reconstitute(peptide, TSrotamer);
                    newPeptide = Rotamer.reconstitute(newPeptide, histidineRotamer);
                    newPeptide = HydrogenBondMutator.mutate(newPeptide);
                    interestingPairPeptides.add(newPeptide);
                }
            //System.out.printf("%d interesting pair peptides\n", interestingPairPeptides.size()); 
            //List<Peptide> minimizedSheets = BetaSheetGenerator.minimizeSheetsInSerial(interestingPairPeptides, 2000, Forcefield.AMOEBA);
            //System.out.printf("%d minimized interesting pair peptides\n", minimizedSheets.size()); 
            List<Peptide> results = addArginine(interestingPairPeptides);
            //List<Peptide> results = BetaSheetGenerator.minimizeSheetsInSerial(argininePeptides, 2000, Forcefield.AMOEBA);
            //System.out.printf("%d arg peptides peptides\n", minimizedSheets.size()); 
            //minimizedSheets = BetaSheetGenerator.minimizeSheetsInSerial(argininePeptides, 2000, Forcefield.AMOEBA);
            //System.out.printf("%d minimized arg peptides\n", minimizedSheets.size()); 
            /*
            // check for proper arg-TS contact
            List<Peptide> resultSheets = new ArrayList<>();
            for (Peptide peptide : minimizedSheets)
                {
                    Atom transitionStateOxygenAtom = null;
                    for (int i=0; i < peptide.sequence.size(); i++)
                        {
                            Residue residue = peptide.sequence.get(i);
                            if ( residue.aminoAcid == AminoAcid.TS )
                                {
                                    for (Atom a : residue.atoms)
                                        {
                                            if ( a.type1 == 428 )
                                                {
                                                    if ( transitionStateOxygenAtom == null )
                                                        transitionStateOxygenAtom = a;
                                                    else
                                                        throw new IllegalArgumentException("duplicate TS O atoms");
                                                }
                                        }
                                }
                        }
                    if ( transitionStateOxygenAtom == null )
                        throw new IllegalArgumentException("didn't find TS oxygen");

                    boolean interesting = false;
                    for (Atom a : peptide.contents)
                        {
                            if ( a.type1 == 209 || a.type1 == 212 )
                                {
                                    if ( Molecule.getDistance(a, transitionStateOxygenAtom) < 2.50 )
                                        {
                                            interesting = true;
                                            break;
                                        }
                                }
                        }
                    if ( interesting )
                        resultSheets.add(peptide);
                }
            */
            if ( results.size() > 0 )
                System.out.printf("found %d solutions\n", results.size());
            interestingPeptides.addAll(results);
            return null;
        }
    }

    /**
     * Searches through a list of atoms and returns the one
     * Throws an exception if the number of matches is not exactly one.
     * @param types the AMOEBA types of the atom we're looking for
     * @param atoms the atoms to search through
     * @return the matching atom
     */
    public static Atom locateSingleAtom(Set<Integer> types, List<Atom> atoms)
    {
        if ( types == null || atoms == null )
            throw new NullPointerException("nulls not allowed");
        Atom resultAtom = null;
        for (Atom a : atoms)
            {
                if ( types.contains(a.type1) )
                    {
                        if ( resultAtom == null )
                            resultAtom = a;
                        else
                            {
                                List<Integer> allTypes = new ArrayList<>();
                                for ( Atom a2 : atoms )
                                    allTypes.add(a2.type1);
                                System.out.println(allTypes);
                                throw new IllegalArgumentException("multiple matches found for " + types.toString());
                            }
                    }
            }
        if ( resultAtom == null )
            throw new NullPointerException("atom not found");
        return resultAtom;
    }

    /**
     * Searches through a list of atoms and returns any that match the specified types.
     * No exceptions are thrown and duplicates are all returned.
     * @param types the AMOEBA types of the atoms we're looking for
     * @param atoms the atoms to search through
     * @return the matching atoms
     */
    public static List<Atom> locateAtoms(Set<Integer> types, List<Atom> atoms)
    {
        if ( types == null || atoms == null )
            throw new NullPointerException("nulls not allowed");
        List<Atom> resultAtoms = new ArrayList<>();;
        for (Atom a : atoms)
            {
                if ( types.contains(a.type1) )
                    resultAtoms.add(a);
            }
        return resultAtoms;
    }

    /** For testing. */
    public static void main(String[] args)
    {
        DatabaseLoader.go();
        List<Peptide> sheets = BetaSheetGenerator.generateSheets(5, 20, 10000, 0.01);
        Collections.sort(sheets);
        int numberOfPeptides = 200;
        System.out.printf("%d beta sheets generated\n", sheets.size());
        List<Peptide> results = new ArrayList<>(numberOfPeptides);
        for (int i=0; i < Math.min(sheets.size(), numberOfPeptides); i++)
            results.add(sheets.get(i));
        
        List<Peptide> interestingPeptides = Collections.synchronizedList(new ArrayList<Peptide>());
        List<Future<Result>> futures = new ArrayList<>();
        for (Peptide peptide : results)
            {
                InterestingJob job = new InterestingJob(peptide, interestingPeptides);
                Future<Result> f = GeneralThreadService.submit(job);
                futures.add(f);
            }
        GeneralThreadService.waitForFutures(futures);

        System.out.printf("\n%d peptides generated\n", interestingPeptides.size());
        Peptide.writeGJFs(interestingPeptides, "test_peptides/result_", 3, 100);
    }
}
