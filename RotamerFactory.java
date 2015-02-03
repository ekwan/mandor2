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
     * Finds the possible rotamers.  This is for "normal" amino acids (i.e., not transition state, histidine, or hairpin).
     * Assumes that the residue is already of the correct type.  This method must return at least one result.
     * If there are no rotamers possible, a rotamer using the existing atoms in the peptide is returned.
     * @param peptide the peptide to analyze
     * @param residue the residue to generate rotamers for
     * @param includeHN if set to true, the backboneHN of this residue will be included in the atoms list of the resulting rotamers
     * @return a list of rotamers that are possible at the specified position
     */
    public static List<Rotamer> generateRotamers(Peptide peptide, Residue residue, boolean includeHN)
    {
        // check this is an appropriate residue
        if ( !peptide.sequence.contains(residue) )
            throw new IllegalArgumentException("residue not in peptide");

        // setup some variables
        int sequenceIndex = peptide.sequence.indexOf(residue);
        String description = residue.description;
        AminoAcid aminoAcid = residue.aminoAcid;
        if ( aminoAcid.chirality == Chirality.D )
            throw new IllegalArgumentException("can't use this method to generate rotamers for " + description);
        List<Rotamer> returnList = new ArrayList<>(); // result rotamers will be placed here

        // get possible chi angles
        List<List<Double>> possibleChis = null;
        if ( residue.description.indexOf("transition_state") > -1 )
            possibleChis = getTransitionStateChis(peptide, residue, includeHN);
        else if ( residue.description.indexOf("histidine") > -1 )
            throw new IllegalArgumentException("not supported");
        else
            possibleChis = RotamerSummarizer.getPossibleRotamers(residue);
        if ( possibleChis.size() == 0 && residue.aminoAcid.rotamerType != AminoAcid.RotamerType.HAS_NO_ROTAMERS && residue.aminoAcid != AminoAcid.TS )
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

        // if there are no rotamers, just return the sidechain atoms
        if ( aminoAcid.rotamerType == AminoAcid.RotamerType.HAS_NO_ROTAMERS )
            {
                // special case for where the amino acid has no rotamers
                List<Double> requestedAngles = ImmutableList.of();  // no torsion angles
                Rotamer singleRotamer = new Rotamer(allAtoms, sequenceIndex, requestedAngles, description);
                returnList.add(singleRotamer);
                return returnList;
            }

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
            }
        return sidechainAtoms;
    }

    /**
     * Returns the atoms that are in the backbone.  Backbone atoms are defined as those that will never move during
     * rotamer packing.  If a position has no rotamers, all atoms at the position are considered part of the backbone.
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
                if ( k >= 0 &&                                                           // must be inside the sequence
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
        System.out.println(sequenceIndex + " : " + HNindices);

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

//    /**
//     * Finds the chi values for a histidine amino acid, under the constraint that the pi-nitrogen must be near a
//     * transition state hydroxyl.
//     */
//    public static List<List<Double>> getHistidineChis(Peptide peptide, Residue residue,
//                                                      List<List<Rotamer>> rotamerSpace, List<Atom> backboneAtoms)
//    {
//        // make a list of interesting atoms at other positions; i.e., the hydroxyl hydrogens
//        int i = peptide.sequence.indexOf(residue);
//        List<Atom> interestingAtoms = new ArrayList<>();
//        List<Atom> interestingOxygenAtoms = new ArrayList<>(); // make a list of the corresponding hydroxyl oxygens
//        for (int j=0; j < rotamerSpace.size(); j++)
//            {
//                if ( i==j )
//                    continue;
//                List<Rotamer> list = rotamerSpace.get(j);
//                for (Rotamer r : list)
//                    {
//                        if ( r.protoAminoAcid.r.description.indexOf("transition_state") > -1 )
//                            {
//                                interestingAtoms.add(r.getInterestingAtom());
//                                boolean partnerFound = false;
//                                for (Atom a : r.atoms)
//                                    {
//                                        if ( a.tinkerAtomType == 402 )
//                                            {
//                                                interestingOxygenAtoms.add(a);
//                                                partnerFound = true;
//                                                break;
//                                            }
//                                    }
//                                if ( !partnerFound )
//                                    throw new IllegalArgumentException("partner not found");
//                            }
//                    }
//            }
//
//
//        // get the rotation axes
//        ProtoTorsion chi1torsion = residue.chis.get(0);
//        ProtoTorsion chi2torsion = residue.chis.get(1);
//
//        // get the sidechain positions
//        Pair<Atom,Atom> prochiralConnection = residue.prochiralConnection;
//        List<Atom> allAtoms = new ArrayList<>(peptide.getHalfGraph(prochiralConnection.getFirst(), prochiralConnection.getSecond()));
//        allAtoms.add(chi1torsion.atom1);
//        allAtoms.add(chi1torsion.atom2);
//
//        // we don't want to check atoms that are part of the backbone for clashes
//        // but just in case there are numerical issues we refer to them by index instead of reference
//        HashSet<Integer> ignoreAtomIndices = new HashSet<>();
//        ignoreAtomIndices.add(allAtoms.indexOf(chi1torsion.atom1));
//        ignoreAtomIndices.add(allAtoms.indexOf(chi1torsion.atom2));
//        ignoreAtomIndices.add(allAtoms.indexOf(chi1torsion.atom3));
//        
//        // the atoms to rotate for each torsion
//        List<Atom> chi1atoms = new ArrayList<>(peptide.getHalfGraph(chi1torsion.atom2, chi1torsion.atom3));
//        List<Atom> chi2atoms = new ArrayList<>(peptide.getHalfGraph(chi2torsion.atom2, chi2torsion.atom3));
//
//        // the indices to rotate for chi2
//        List<Integer> chi2indices = new ArrayList<>();
//        for (Atom a : chi2atoms)
//            {
//                int atomIndex = allAtoms.indexOf(a);
//                chi2indices.add(atomIndex);
//            }
//
//        // the indices of the torsion for chi2
//        List<Integer> chi2torsionIndices = new ArrayList<>();
//        chi2torsionIndices.add(allAtoms.indexOf(chi2torsion.atom1));
//        chi2torsionIndices.add(allAtoms.indexOf(chi2torsion.atom2));
//        chi2torsionIndices.add(allAtoms.indexOf(chi2torsion.atom3));
//        chi2torsionIndices.add(allAtoms.indexOf(chi2torsion.atom4));
//
//        // get the index of the histidine pi-nitrogen
//        Integer atomNindex = null; 
//        for (Atom a : allAtoms)
//            {
//                if ( a.tinkerAtomType == 126 || a.tinkerAtomType == 135 )
//                    {
//                        atomNindex = allAtoms.indexOf(a);
//                        break;
//                    }
//            }
//        if ( atomNindex == null )
//            throw new NullPointerException("null histidine N index");
//
//        // rotate chi1 and chi2 on a grid
//        // check for the desired contact
//        List<List<Double>> returnList = new ArrayList<>();
//        double stepSize = 360.0 / (TS_GRID_SIZE-1.0);
//        double chi1 = -180.0;
//        int sequenceIndex = peptide.sequence.indexOf(residue);
//        if ( sequenceIndex == -1 )
//            throw new IllegalArgumentException("bad sequence index");
//        for (int j=0; j < TS_GRID_SIZE; j++)
//            {
//                List<Atom> thisAtoms = setDihedral(allAtoms, chi1atoms, chi1torsion, chi1);
//                chi2atoms.clear();
//                for (Integer index : chi2indices)
//                    chi2atoms.add(thisAtoms.get(index));
//                Atom atom1 = thisAtoms.get(chi2torsionIndices.get(0));
//                Atom atom2 = thisAtoms.get(chi2torsionIndices.get(1));
//                Atom atom3 = thisAtoms.get(chi2torsionIndices.get(2));
//                Atom atom4 = thisAtoms.get(chi2torsionIndices.get(3));
//                chi2torsion = new ProtoTorsion(atom1, atom2, atom3, atom4);
//
//                double chi2 = -180.0;
//                for (int k=0; k < TS_GRID_SIZE; k++)
//                    {
//                        // try this rotamer
//                        List<Atom> thisAtoms2 = setDihedral(thisAtoms, chi2atoms, chi2torsion, chi2);
//
//                        // get the oxygen
//                        Atom atomN = thisAtoms2.get(atomNindex);
//
//                        // check that this could become an interesting pair
//                        boolean acceptable = false;
//
//                        for ( Atom a : interestingAtoms )
//                            {
//                                double distance = Vector3D.distance(atomN.position, a.position);
//                                if ( distance > RotamerPacker.HBONDED_MIN_DISTANCE && distance < RotamerPacker.INTERESTING_DISTANCE )
//                                    {
//                                        Integer index = interestingAtoms.indexOf(a); 
//                                        Atom oxygenAtom = interestingOxygenAtoms.get(index);
//                                        double angle = Molecule.getAngle(oxygenAtom, a, atomN); 
//                                        if ( angle > HBONDED_MINIMUM_ANGLE )
//                                            {
//                                                acceptable = true;
//                                                break;
//                                            }
//                                    }
//                            }
//
//                        if ( acceptable )
//                            {
//                                // check for backbone atom clashes
//                                HashSet<Atom> ignoreAtoms = new HashSet<>();
//                                for (Integer index : ignoreAtomIndices)
//                                    ignoreAtoms.add(thisAtoms2.get(index));
//
//                                clashCheck:
//                                for ( Atom a1 : thisAtoms2 )
//                                    {
//                                        if ( ignoreAtoms.contains(a1) )
//                                            continue;
//                                        for (Atom a2 : backboneAtoms )
//                                            {
//                                                double distance = Vector3D.distance(a1.position, a2.position);
//                                                if ( distance < RotamerPacker.HBONDED_MIN_DISTANCE )
//                                                    {
//                                                        acceptable = false;
//                                                        break clashCheck;
//                                                    }
//                                            }
//                                    }
//                        
//                                // if acceptable, add it to the list of chis to return
//                                if ( acceptable )
//                                    returnList.add(ImmutableList.of(chi1,chi2));
//                            }
//
//                        // increment chi2
//                        chi2 += stepSize;
//                    }
//                // increment chi1
//                chi1 += stepSize;
//            }
//        //System.out.println("his size: " + returnList.size());
//        //System.out.println("his: " + returnList);
//        return returnList;
//    }
//
//    /**
//     * Returns the interesting atom in this Rotamer.
//     * Only the first hit is returned.  Intended to be used on structures that are not using the close
//     * contact atom types.  See atom typing.cdx.
//     * @return the interesting atom if present, null otherwise
//     */
//    public Atom getInterestingAtom()
//    {
//        for (Atom a : atoms)
//            {
//                if ( a.tinkerAtomType == 126 || a.tinkerAtomType == 130 ||   // histidine pi nitrogens
//                     a.tinkerAtomType == 401                               ) // methanol HO
//                    return a;
//            }
//        return null;
//    }

    /** For testing. */
    public static void main(String[] args)
    {
        DatabaseLoader.go();
        List<Peptide> sheets = BetaSheetGenerator.generateSheets(5, 10, 10000, 0.01);
        Collections.sort(sheets);
        Peptide peptide = sheets.get(0);
        
        int TSindex = 4; // put the TS at this sequence index
        ProtoAminoAcid tsTemplate = ProtoAminoAcidDatabase.getTemplate("ts1");
        peptide = SidechainMutator.mutateSidechain(peptide, peptide.sequence.get(TSindex), tsTemplate);
        new GaussianInputFile(peptide).write("test_peptides/peptide.gjf");

        List<Rotamer> rotamers = generateRotamers(peptide, peptide.sequence.get(TSindex), true);
        for (int i=0; i < rotamers.size(); i++)
            {
                Rotamer rotamer = rotamers.get(i);
                Peptide newPeptide = Rotamer.reconstitute(peptide, rotamer);
                GaussianInputFile f = new GaussianInputFile(newPeptide);
                String filename = String.format("test_peptides/rotamer_%02d.gjf", i);
                f.write(filename);
                System.out.printf("%02d %s\n", i, rotamer.toString());
            }
    }
}
