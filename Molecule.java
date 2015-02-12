import java.util.*;
import java.io.*;
import org.apache.commons.math3.geometry.euclidean.threed.*;
import org.jgrapht.*;
import org.jgrapht.graph.*;
import org.jgrapht.alg.*;
import com.google.common.collect.*;
import Jama.*;

/**
 * Represents a molecule.  This class is effectively immutable.
 */
public class Molecule implements Immutable, Serializable
{
    /** For serialization. */
    public static final long serialVersionUID = 1L;

    /** Name of the molecule. */
    public final String name;

    /** The atoms in the molecule. */
    public final List<Atom> contents;

    /** The connectivity graph. */
    public final SimpleWeightedGraph<Atom, DefaultWeightedEdge> connectivity;

    /**
     * Factory method to create a molecule given a map of old atoms to new atoms.  Should be used
     * to move atoms.
     * @param atomMap a map from old atoms to new atoms (does not have to include all atoms)
     */
    public Molecule moveAtoms(Map<Atom,Atom> atomMap)
    {
        // copy the list of vertices
        List<Atom> newContents = new ArrayList<Atom>(contents.size());
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
                if (! newContents.contains(fromAtom) || ! newContents.contains(toAtom))
                    throw new IllegalArgumentException("edge not in graph   FromAtom: " + getAtomString(fromAtom) +  " ToAtom: " + getAtomString(toAtom));
                DefaultWeightedEdge newEdge = newConnectivity.addEdge(fromAtom,toAtom);
                newConnectivity.setEdgeWeight(newEdge, bondOrder);
            }

        // return result
        return new Molecule(name, newContents, newConnectivity);
    }

    /** Adds the specified shift to all the atoms. */
    public Molecule shift(Vector3D shift)
    {
        Map<Atom,Atom> atomMap = new HashMap<>();
        for (Atom a : contents)
            atomMap.put(a, a.moveAtom(a.position.add(shift)));
        return moveAtoms(atomMap);
    }

    /**
     * A shallow clone of a connectivity graph.
     * @param oldConnectivity the old connectivity graph
     * @return the clone
     */
    public static SimpleWeightedGraph<Atom,DefaultWeightedEdge> cloneConnectivity(SimpleWeightedGraph<Atom,DefaultWeightedEdge> oldConnectivity)
    {
        // populate a new connectivity graph
        SimpleWeightedGraph<Atom,DefaultWeightedEdge> newConnectivity = new SimpleWeightedGraph<Atom,DefaultWeightedEdge>(DefaultWeightedEdge.class);
        for (Atom oldAtom : oldConnectivity.vertexSet())
            newConnectivity.addVertex(oldAtom);
        for (DefaultWeightedEdge e : oldConnectivity.edgeSet())
            {
                Double bondOrder = oldConnectivity.getEdgeWeight(e);
                Atom fromAtom    = oldConnectivity.getEdgeSource(e);
                Atom toAtom      = oldConnectivity.getEdgeTarget(e);
                DefaultWeightedEdge newEdge = newConnectivity.addEdge(fromAtom,toAtom);
                newConnectivity.setEdgeWeight(newEdge, bondOrder);
            }
        return newConnectivity;
    }

    /**
     * Factory method to create a new molecule by shifting this one to have its center at the origin.
     * @return this minus the barycenter of this 
     */
    public Molecule normalize()
    {
        Vector3D barycenter = Vector3D.ZERO;
        for (Atom a : contents) barycenter = barycenter.add(a.position);
        barycenter = barycenter.scalarMultiply(1.0/contents.size());
        return shift(barycenter.negate());
    }

    /**
     * Returns the centroid of this molecule.
     * @return the centroid vector
     */
    public Vector3D getCentroid()
    {
        Vector3D barycenter = Vector3D.ZERO;
        for (Atom a : contents) barycenter = barycenter.add(a.position);
        barycenter = barycenter.scalarMultiply(1.0/contents.size());
        return barycenter;
    }

    /**
     * Constructs a new molecule.  The connectivity graph will be used directly (i.e., no defensive copy is made).
     * @param name the name of the molecule
     * @param contents the atoms in the molecule
     * @param connectivity the connectivity graph
     */
    public Molecule(String name, List<Atom> contents, SimpleWeightedGraph<Atom,DefaultWeightedEdge> connectivity)
    {
	    this.name = name;
	    this.contents = ImmutableList.copyOf(contents);
        this.connectivity = connectivity;
    }

    /**
     * Returns the number of the atom.  Answers are given in 1,2,...n where n is
     * the number of atoms in the molecule.
     * @param atom an Atom in this Molecule
     * @return the requested atom number (0 if it is not in the molecule)
     */
    public int getAtomNumber(Atom atom)
    {
        return contents.indexOf(atom) + 1;
    }

    /**
     * Returns the element and atom number of an atom.  e.g., C5.
     * Throws an IllegalArgumentException if the Atom is not in the molecule.
     * @param atom the atom whose string is desired
     * @return the description of the atom
     */
    public String getAtomString(Atom atom)
    {
        if ( atom == null)
            throw new NullPointerException("atom cannot be null");

        int atomNumber = getAtomNumber(atom);
        
        if ( atomNumber == 0 )
            throw new IllegalArgumentException("cannot print out atom because it is not in the molecule!\n" + atom.toString());
        return atom.element.symbol + atomNumber;
    }

    /** For debugging. */
    public String getClosestAtomString(Atom atom)
    {
        double closestDistance = 100.0;
        Atom closestAtom = null;
        for (Atom a : contents)
            {
                if ( a.type1 != atom.type1 )
                    continue;
                double distance = Vector3D.distance(atom.position, a.position);
                if ( distance < closestDistance )
                    {
                        closestDistance = distance;
                        closestAtom = a;
                    }
            }
        if ( closestDistance > 0.1 || closestAtom == null )
            return atom.toString();
        else
            return getAtomString(closestAtom);
    }

    /**
     * Given a single atom, returns a set of all connected atoms, including the
     * specified atom.  Uses a breadth-first search algorithm.
     * @param startingAtom the atom to start exploring the connectivity from
     * @return the atoms in the subgraph of the startingAtom
     */
    public Set<Atom> exploreGraph(Atom startingAtom)
    {
        if ( ! contents.contains(startingAtom) )
            throw new IllegalArgumentException("cannot search the connectivity graph because the specified atom is not in this molecule");
        Set<Atom> returnSet = new HashSet<>();
        LinkedList<Atom> searchQueue = new LinkedList<>();
        searchQueue.add(startingAtom);

        // breadth-first search
        while (searchQueue.size() > 0)
            {
                Atom currentNode = searchQueue.remove();
                for (Atom a : getAdjacentAtoms(currentNode))
                    {
                        if ( ! returnSet.contains(a) )
                            {
                                returnSet.add(a);
                                searchQueue.add(a);
                            }
                    }
            }
        return returnSet;
    }

    /**
     * Given a bond between includeAtom and excludeAtom, returns a set of atoms
     * containing all the atoms on the includeAtom side of the bond, including
     * includeAtom.  Will throw an exception if includeAtom and excludeAtom are not
     * directly bonded.  Will throw an exception if includeAtom and excludeAtom
     * form a ring.
     * @param excludeAtom this atom will not be included in the result
     * @param includeAtom this atom will be included in the result
     * @return the atoms on the includeAtom side of the graph
     */
    public Set<Atom> getHalfGraph(Atom excludeAtom, Atom includeAtom)
    {
        // if these atoms are not directly bonded, then return an empty set
        if ( !directlyConnected(excludeAtom, includeAtom) )
            throw new IllegalArgumentException(String.format("atoms %s and %s not directly connected", getAtomString(excludeAtom), getAtomString(includeAtom))); 

        // preform a breadth-first search of one branch of the graph only
        Set<Atom> returnSet = new HashSet<Atom>();
        LinkedList<Atom> searchQueue = new LinkedList<Atom>();
        Set<Atom> searched = new HashSet<>();
        for (Atom a : getAdjacentAtoms(includeAtom))
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
                Set<Atom> adjacent = getAdjacentAtoms(currentNode);
                adjacent.remove(lastNode);
                for (Atom a : adjacent)
                    {
                        if ( a == excludeAtom)
                            {
                                // if the excluded atom is found, this is a ring!
                                GaussianInputFile gjf = new GaussianInputFile(this);
                                gjf.write("error.gjf");
                                throw new IllegalArgumentException("includeAtom " + getAtomString(includeAtom) +
                                                                   " and excludeAtom " + getAtomString(excludeAtom) + " cannot form a ring!");
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
     * {@link #getHalfGraph(Atom,Atom)}
     * @param index1 the index of the excludeAtom
     * @param index2 the index of the includeAtom
     * @return the atoms on the includeAtom side of the graph
     */
    public Set<Atom> getHalfGraph(int index1, int index2)
    {
        Atom atom1 = contents.get(index1);
        Atom atom2 = contents.get(index2);
        if ( atom1 == null || atom2 == null )
            throw new NullPointerException("atoms not found in this molecule for half graph call");
        return getHalfGraph(atom1, atom2);
    }
    
    /**
     * Returns atom indices instead of atoms.
     * {@link #getHalfGraph(Atom,Atom)}
     * @param excludeAtom this atom will not be included in the result
     * @param includeAtom this atom will be included in the result
     * @return the atom indices on the includeAtom side of the graph
     */
    public Set<Integer> getHalfGraphIndices(Atom excludeAtom, Atom includeAtom)
    {
        Set<Atom> set = getHalfGraph(excludeAtom, includeAtom);
        Set<Integer> returnSet = new HashSet<>();
        for (Atom a : set)
            returnSet.add(contents.indexOf(a));
        return returnSet;
    }

    /**
     * Returns atom indices instead of atoms.
     * {@link #getHalfGraph(Atom,Atom)}
     * @param index1 the excludeAtom index
     * @param index2 the includeAtom index
     * @return the atom indices on the includeAtom side of the graph
     */

    public Set<Integer> getHalfGraphIndices(int index1, int index2)
    {
        return getHalfGraphIndices(contents.get(index1), contents.get(index2));
    }

    /**
     * Determines whether atom1 and atom2 share an edge (i.e., are bonded).
     * No exception is thrown if these atoms aren't in the graph.
     * @param atom1 test whether this atom is connected to the other atom 
     * @param atom2 test whether this atom is connected to the other atom
     * @return true if the atoms are bonded
     */
    public boolean directlyConnected(Atom atom1, Atom atom2)
    {
        DefaultWeightedEdge e = connectivity.getEdge(atom1,atom2);
        if ( e == null )
            return false;
        return true;
    }

    /**
     * {@link #directlyConnected(Atom,Atom)}
     * @param i index of atom1
     * @param j index of atom2
     * @return true if atoms are bonded
     */
    public boolean directlyConnected(int i, int j)
    {
        return directlyConnected(contents.get(i), contents.get(j));
    }

    /**
     * Returns the distance between two atoms.
     * @param atom1 one of the two atoms
     * @param atom2 one of the two atoms
     * @return the distance between the atoms in angstroms
     */
    public static double getDistance(Atom atom1, Atom atom2)
    {
        return Vector3D.distance(atom1.position, atom2.position);
    }

    /**
     * Returns the angle between three atoms.
     * @param atom1 one of the three atoms
     * @param atom2 one of the three atoms
     * @param atom3 one of the three atoms
     * @return the angle between the atoms in degrees
     */
    public static double getAngle(Atom atom1, Atom atom2, Atom atom3)
    {
        return getAngle(atom1.position, atom2.position, atom3.position);
    }

    /**
     * Returns the v1-v2-v3 angle.
     * @param v1 the first point
     * @param v2 the second point
     * @param v3 the third point
     * @return the angle in degrees
     */
    public static double getAngle(Vector3D v1, Vector3D v2, Vector3D v3)
    {
        Vector3D v1prime = v1.subtract(v2);
        Vector3D v3prime = v3.subtract(v2);
        return Math.toDegrees(Vector3D.angle(v1prime, v3prime));
    }

    /**
     * Returns the angle between three atoms.
     * @param index1 the index of atom1
     * @param index2 the index of atom2
     * @param index3 the index of atom3
     * @return the angle in degrees
     */
    public double getAngle(int index1, int index2, int index3)
    {
        return getAngle(contents.get(index1), contents.get(index2), contents.get(index3));
    }

    /**
     * Moves the group associated with atom2 to the specified distance.
     * Motion occurs along the atom1-atom2 bond vector.  Note that this returns a new
     * molecule.  No checks are made.
     * @param atom1 this atom will be held fixed
     * @param atom2 this atom and anything connected to it will be moved
     * @param requestedDistance the requested distance in Angstroms
     * @return a new Molecule containing the same connectivity but new positions
     */
    public Molecule setDistance(Atom atom1, Atom atom2, double requestedDistance)
    {
        // determine which atoms have to be moved
        Set<Atom> toBeMoved = getHalfGraph(atom1, atom2);

        // determine how much to move the atoms
        Vector3D oldPosition1 = atom1.position;
        Vector3D oldPosition2 = atom2.position;
        double currentDistance = getDistance(atom1, atom2);

        Vector3D translateVector = oldPosition2.subtract(oldPosition1);
        double scaling = (requestedDistance - currentDistance)/currentDistance;

        Vector3D requiredTranslation = translateVector.scalarMultiply(scaling);

        Map<Atom,Atom> atomMap = new HashMap<>();
        for (Atom oldAtom : toBeMoved)
            {
                Vector3D oldPosition = oldAtom.position;
                Vector3D newPosition = oldPosition.add(requiredTranslation);
                Atom newAtom         = oldAtom.moveAtom(newPosition);
                atomMap.put(oldAtom, newAtom);
            }
       return moveAtoms(atomMap); 
    }

    /**
     * Moves the group associated with atom2 to the specified distance.
     * Motion occurs along the atom1-atom2 bond vector.  Note that this returns a new
     * molecule.  No checks are made.  Indices are 0, 1, ..., n-1 (the as the contents field).
     * @param index1 this atom index will be held fixed
     * @param index2 this atom indexand anything connected to it will be moved
     * @param requestedDistance the requested distance in Angstroms
     * @return a new Molecule containing the same connectivity but new positions
     */
    public Molecule setDistance(int index1, int index2, double requestedDistance)
    {
        return setDistance(contents.get(index1), contents.get(index2), requestedDistance);
    }

    /**
     * Rotates the atom1-atom2-atom3 angle, moving only atom3 and anything in its
     * attached subgraph.  No checks.
     * @param atom1 will not be moved
     * @param atom2 will not be moved
     * @param atom3 will be moved
     * @param theta rotation in degrees
     * @return a new Molecule containing the same connectivity but new positions
     */
    public Molecule rotateAngle(Atom atom1, Atom atom2, Atom atom3, double theta)
    {
        // figure out which atoms to move
        Set<Atom> toBeMoved = getHalfGraph(atom2, atom3);

        // create atom map
        LinkedHashMap<Atom,Atom> atomMap = new LinkedHashMap<Atom,Atom>();

        // move everything to put atom2 at the origin
        Vector3D v1 = null;
        Vector3D v3 = null;
        for (Atom a : contents)
            {
                Vector3D oldPosition = a.position;
                Vector3D newPosition = oldPosition.subtract(atom2.position);
                if ( toBeMoved.contains(a) )
                    atomMap.put(a,a.moveAtom(newPosition));
                if ( a == atom1 )
                    v1 = newPosition;
                else if ( a == atom3 )
                    v3 = newPosition;
            }
        
        // form the rotation axis and matrix
        Vector3D rotationAxis = Vector3D.crossProduct(v1, v3);
        Rotation rotation = new Rotation(rotationAxis, Math.toRadians(theta));

        // apply rotation and undo translation
        LinkedHashMap<Atom,Atom> atomMap2 = new LinkedHashMap<Atom,Atom>();
        for (Atom a : atomMap.keySet())
            {
                Vector3D oldPosition = atomMap.get(a).position;
                Vector3D newPosition = rotation.applyTo(oldPosition);
                newPosition = newPosition.add(atom2.position);
                atomMap2.put(a,a.moveAtom(newPosition));
            }

        // create new Molecule
        return moveAtoms(atomMap2);
    }

    /**
     * Method alias.  Indices are 0, 1, ..., n-1.  No checks.
     */
    public Molecule rotateAngle(int i, int j, int k, double theta)
    {
        return rotateAngle(contents.get(i), contents.get(j), contents.get(k), theta);
    }

    /**
     * Set the atom1-atom2-atom3 angle to theta degrees, moving atom3 and its subgraph only.
     * New molecule returned.  No checks.
     * @param atom1 not moved
     * @param atom2 not moved
     * @param atom3 moved
     * @param theta desired angle in degrees
     */
    public Molecule setAngle(Atom atom1, Atom atom2, Atom atom3, double theta)
    {
        double currentAngle = getAngle(atom1, atom2, atom3);
        double requiredRotation = theta - currentAngle;
        return rotateAngle(atom1, atom2, atom3, requiredRotation);
    }

    /**
     * Method alias.  Indices are 1,2,...,n.  No checks.
     */
    public Molecule setAngle(int i, int j, int k, double theta)
    {
        return setAngle(contents.get(i), contents.get(j), contents.get(k), theta);
    }

    /**
     * Returns a new Molecule with a rotated dihedral.
     * Note that the old ProtoTorsion will no longer point to the new Molecule.
     * @param protoTorsion the torsion to rotate
     * @param theta the desired dihedral angle in degrees
     * @return the new molecule
     */
    public Molecule setDihedral(ProtoTorsion protoTorsion, double theta)
    {
        // get fields
        Atom atom1 = protoTorsion.atom1;
        Atom atom2 = protoTorsion.atom2;
        Atom atom3 = protoTorsion.atom3;
        Atom atom4 = protoTorsion.atom4;
        Set<Atom> atomsToRotate = getHalfGraph(atom2,atom3);

        // determine how much rotation is needed
        double currentDihedralAngle = protoTorsion.getDihedralAngle();
        double requiredRotation = currentDihedralAngle - theta;
        if ( Math.abs(requiredRotation) < 0.001 )
            return this;

        Map<Atom,Atom> atomMap = new HashMap<>();
        
        // move atom 3 to the origin
        // define the rotation axis as the vector from atom3 (now at origin) to atom2
        Vector3D rotationAxis = null;
        for (Atom a : contents)
            {
                Vector3D oldPosition = a.position;
                Vector3D newPosition = oldPosition.subtract(atom3.position);
                if ( atomsToRotate.contains(a) )
                    atomMap.put(a,a.moveAtom(newPosition));
                if ( a == atom2 )
                    rotationAxis = newPosition;
            }

        // rotate the atoms and make a new atom map
        Map<Atom,Atom> atomMap2 = new HashMap<>();
        Rotation rotation = new Rotation(rotationAxis, Math.toRadians(requiredRotation));
        for (Atom a : atomMap.keySet())
            {
                // update rotation
                Vector3D oldPosition = atomMap.get(a).position;
                Vector3D newPosition = rotation.applyTo(oldPosition);
                
                // undo translation
                newPosition = newPosition.add(atom3.position);

                // update map
                atomMap2.put(a, a.moveAtom(newPosition));
            }

        // return new Molecule
        return moveAtoms(atomMap2);
    }

    /**
     * Returns a new Molecule with a rotated dihedral.  Alias method.
     * Note that the old IndexTorsion will still be valid for the new Molecule.
     * @param indexTorsion the torsion to rotate
     * @param theta the angle to set the torsion to in degrees
     */
    public Molecule setDihedral(IndexTorsion indexTorsion, double theta)
    {
        // figure out which atoms to move
        Set<Atom> atomsToRotate = new HashSet<>();
        for (Integer i : indexTorsion.atomIndicesToRotate)
            atomsToRotate.add(contents.get(i));

        // get prototorsion
        ProtoTorsion protoTorsion = indexTorsion.getProtoTorsion(this);
        Atom atom1 = protoTorsion.atom1;
        Atom atom2 = protoTorsion.atom2;
        Atom atom3 = protoTorsion.atom3;
        Atom atom4 = protoTorsion.atom4;
        
        // determine how much rotation is needed
        double currentDihedralAngle = protoTorsion.getDihedralAngle();
        double requiredRotation = currentDihedralAngle - theta;
        if ( Math.abs(requiredRotation) < 0.001 )
            return this;

        Map<Atom,Atom> atomMap = new HashMap<>();
        
        // move atom 3 to the origin
        // define the rotation axis as the vector from atom3 (now at origin) to atom2
        Vector3D rotationAxis = null;
        for (Atom a : contents)
            {
                Vector3D oldPosition = a.position;
                Vector3D newPosition = oldPosition.subtract(atom3.position);
                if ( atomsToRotate.contains(a) )
                    atomMap.put(a,a.moveAtom(newPosition));
                if ( a == atom2 )
                    rotationAxis = newPosition;
            }

        // rotate the atoms and make a new atom map
        Map<Atom,Atom> atomMap2 = new HashMap<>();
        Rotation rotation = new Rotation(rotationAxis, Math.toRadians(requiredRotation));
        for (Atom a : atomMap.keySet())
            {
                // update rotation
                Vector3D oldPosition = atomMap.get(a).position;
                Vector3D newPosition = rotation.applyTo(oldPosition);
                
                // undo translation
                newPosition = newPosition.add(atom3.position);

                // update map
                atomMap2.put(a, a.moveAtom(newPosition));
            }

        // return new Molecule
        return moveAtoms(atomMap2);
    }

    /**
     * Creates a new Molecule where atom2 and its subgraph have been moved to
     * make atom1 sp2-hybridized (bond angles set at 120 degrees).
     * Note that the new center will be sp2, but could have distorted torsion angles.
     * @param atom1 the atom to be adjusted to sp2
     * @param atom2 the group to be moved
     * @param forceAngle true if we want to force the atom1alpha-atom1-atom1beta angle to 120 (safe for non-prolines)
     * @return a new Molecule with adjusted hybridization
     */
    public Molecule set_sp2(Atom atom1, Atom atom2, boolean forceAngle)
    {
        // note current bond length
        double currentLength = getDistance(atom1, atom2);

        // get the neighbors of atom1
        // call them atom1alpha and atom1beta
        List<Atom> atom1neighbors = new LinkedList<>(getAdjacentAtoms(atom1));
        if ( atom1neighbors.size() != 3 )
            throw new IllegalArgumentException("expected 3 neighbors for atom 1, found " + atom1neighbors.size());
        else if ( ! atom1neighbors.contains(atom2) )
            throw new IllegalArgumentException("atoms are not adjacent");
        atom1neighbors.remove(atom2);
        Atom atom1alpha = atom1neighbors.get(0);
        Atom atom1beta = atom1neighbors.get(1);
        int atom1alphaIndex = contents.indexOf(atom1alpha);
        int atom1betaIndex = contents.indexOf(atom1beta);
        int atom1index = contents.indexOf(atom1);
        int atom2index = contents.indexOf(atom2);

        // force the existing bond angle to 120 degrees
        Molecule newMolecule = this;
        if (forceAngle)
            newMolecule = setAngle(atom1alpha, atom1, atom1beta, 120.0);
        //System.out.println(getAtomString(atom1alpha));
        //System.out.println(getAtomString(atom1));
        //System.out.println(getAtomString(atom1beta));

        // translate atom1 to the origin
        List<Vector3D> newPositions = new LinkedList<>();
        for (Atom a : newMolecule.contents)
            {
                Vector3D oldPosition = a.position;
                Vector3D newPosition = oldPosition.subtract(atom1.position);
                newPositions.add(newPosition);
            }

        // get unit vectors for atom1alpha and atom1beta
        Vector3D atom1alphaPosition = newPositions.get(atom1alphaIndex).normalize();
        Vector3D atom1betaPosition = newPositions.get(atom1betaIndex).normalize();

        // get the cross product of atom1alpha and atom1beta
        Vector3D atom1crossPosition = Vector3D.crossProduct(atom1alphaPosition, atom1betaPosition).normalize();

        // find the linear transformation matrix that rotates atom1alphaPosition to a=(-sqrt(3)/2, 0.5, 0.0)
        // and atom1betaPosition to b=(-0.5, sqrt(3)/2, 0.0). If these are two points of an equilateral triangle,
        // the third is at c=(1,0,0).

        // solve the simultaneous matrix equations:
        // T atom1alphaPosition = a
        // T atom1betaPosition = b
        // T atom1cross = c
        //
        // call atom1alphaPosition A, atom1betaPosition B, and atom1cross C.
        // this is equivalent to solving:
        //
        // T [ ABC ] = [ abc ]
        //
        // where this means concatenated column vectors. Therefore,
        // T = [ abc ] [ ABC ]^-1

        double[][] preMatrix_abc = { { -0.5, Math.sqrt(3.0)/2.0, 0.0 }, {-0.5, -1.0 * Math.sqrt(3.0)/2.0, 0.0}, {0.0, 0.0, 1.0} };
        Matrix matrix_abc = new Matrix(preMatrix_abc);
        matrix_abc = matrix_abc.transpose();

        double[][] preMatrix_ABC = { { atom1alphaPosition.getX(), atom1alphaPosition.getY(), atom1alphaPosition.getZ() },
                                     { atom1betaPosition.getX(), atom1betaPosition.getY(), atom1betaPosition.getZ() },
                                     { atom1crossPosition.getX(), atom1crossPosition.getY(), atom1crossPosition.getZ() } };
        Matrix matrix_ABC = new Matrix(preMatrix_ABC);
        matrix_ABC = matrix_ABC.transpose();

        Matrix matrix_ABC_inverse = matrix_ABC.inverse();
        Matrix T = matrix_abc.times(matrix_ABC_inverse);
        Matrix Tinverse = T.inverse();

        // apply the inverse of T to (1,0,0) to get the third vertex of the triangle
        double[][] preMatrix_c = { { 1.0, 0.0, 0.0 } };
        Matrix matrix_c = new Matrix(preMatrix_c);
        matrix_c = matrix_c.transpose();

        Matrix thirdVertex = Tinverse.times(matrix_c);
        Vector3D thirdVertexPosition = new Vector3D( thirdVertex.get(0, 0), thirdVertex.get(1, 0), thirdVertex.get(2, 0) );

        //double angle1 = getAngle(newPositions.get(atom1alphaNumber-1), newPositions.get(atom1number-1), newPositions.get(atom1betaNumber-1));
        //double angle2 = getAngle(newPositions.get(atom1betaNumber-1), newPositions.get(atom1number-1), newPositions.get(atom2number-1));
        //double angle3 = getAngle(newPositions.get(atom2number-1), newPositions.get(atom1number-1), newPositions.get(atom1alphaNumber-1));
        //System.out.println(angle1);
        //System.out.println(angle2);
        //System.out.println(angle3);

        // calculate the necessary rotation
        Vector3D atom2position = newPositions.get(atom2index);
        Vector3D rotationAxis = Vector3D.crossProduct( thirdVertexPosition, atom2position );
        double requiredTheta = Vector3D.angle( thirdVertexPosition, atom2position );
        Rotation fixRotation = new Rotation( rotationAxis, -1.0 * requiredTheta );

        // determine which atoms should be moved
        Set<Integer> atomIndicesToMove = getHalfGraphIndices(atom1, atom2);
        List<Vector3D> newPositions2 = new LinkedList<>();

        for (int i=0; i < newPositions.size(); i++)
            {
                if ( atomIndicesToMove.contains(i) )
                    {
                        // rotate this atom
                        Vector3D oldPosition = newPositions.get(i);
                        Vector3D newPosition = fixRotation.applyTo(oldPosition);
                        newPositions2.add(newPosition);
                    }
                else
                    {
                        // do not rotate this atom
                        newPositions2.add( newPositions.get(i) );
                    }
            }

        // undo translation
        List<Vector3D> newPositions3 = new LinkedList<>();
        for (Vector3D v : newPositions2)
            newPositions3.add(v.add(atom1.position));

        /*double angle1 = getAngle(newPositions3.get(atom1alphaIndex), newPositions3.get(atom1index), newPositions3.get(atom1betaIndex));
        double angle2 = getAngle(newPositions3.get(atom1betaIndex), newPositions3.get(atom1index), newPositions3.get(atom2index));
        double angle3 = getAngle(newPositions3.get(atom2index), newPositions3.get(atom1index), newPositions3.get(atom1alphaIndex));
        System.out.println(angle1);
        System.out.println(angle2);
        System.out.println(angle3);*/
        //System.out.println(atom1alphaNumber + " " + atom1number + " " + atom1betaNumber);
        //System.out.println(atom1betaNumber + " " + atom1number + " " + atom2number);
        //System.out.println(atom2number + " " + atom1number + " " + atom1alphaNumber);

        // create new atom map
        Map<Atom,Atom> newAtomMap = new HashMap<>();
        for (int i=0; i < contents.size(); i++)
            {
                Atom oldAtom = contents.get(i);
                Atom newAtom = oldAtom.moveAtom( newPositions3.get(i) );
                if ( !oldAtom.equals(newAtom) )
                    newAtomMap.put(oldAtom, newAtom);
            }

        Molecule rotatedMolecule = moveAtoms(newAtomMap);

        
        //for (Atom key : newAtomMap.keySet())
        //    System.out.println(String.format("%s %s  :  %s %s", getAtomString(key), key, rotatedMolecule.getAtomString(newAtomMap.get(key)), newAtomMap.get(key)));
        
        // set bond length
        Molecule returnMolecule = rotatedMolecule.setDistance(atom1index, atom2index, currentLength);
        //System.out.println(returnMolecule.getAngle(atom1alphaNumber,atom1number,atom1betaNumber));
        //System.out.println(returnMolecule.getAngle(atom1betaNumber,atom1number,atom2number));
        //System.out.println(returnMolecule.getAngle(atom2number,atom1number,atom1alphaNumber));
        
        return returnMolecule;
    }

    public Molecule set_sp2(Atom atom1, Atom atom2)
    {
        return set_sp2(atom1, atom2, true);
    }

    /**
     * Returns the bonded neighbors of includeAtom.  Does not include includeAtom itself.
     * @param includeAtom the atom whose neighbors are to be searched
     * @return the Atoms adjacent to includeAtom
     */
    public Set<Atom> getAdjacentAtoms(Atom includeAtom)
    {
        Set<Atom> returnSet = new HashSet<Atom>();
        if ( ! connectivity.containsVertex(includeAtom) )
            throw new IllegalArgumentException("includeAtom must be within this connectivity graph!");
        for (DefaultWeightedEdge e : connectivity.edgesOf(includeAtom))
            {
                returnSet.add(connectivity.getEdgeSource(e));
                returnSet.add(connectivity.getEdgeTarget(e));
            }
        returnSet.remove(includeAtom);
        return(returnSet);
    }

    /**
     * Checks if atom1 and atom2 are more than two bonds apart.
     * @param atom1 the first atom
     * @param atom2 the second atom
     * @return true if atom1 and atom2 are separated by three or more bonds
     */
    public boolean areSeparated(Atom atom1, Atom atom2)
    {
        Set<Atom> atom1neighbors = getAdjacentAtoms(atom1);
        
        // check if direct neighbors
        if (atom1neighbors.contains(atom2))
            return false;
        Set<Atom> atom2neighbors = getAdjacentAtoms(atom2);
       
        // check if geminal
        for (Atom a : atom1neighbors)
            if (atom2neighbors.contains(a))
                return false;
        return true;
    }
 
    /**
     * Checks if atoms are too close in a molecule, given another molecule
     * whose atoms we know are not too close.  Intended for assessing the result
     * of dihedral changes.  Does not assume the molecules have the same composition.
     * @param oldMolecule the molecule this molecule was modified from
     * @return true if there is at least one atom that is too close to another atom
     */
    public boolean tooClose(Molecule oldMolecule)
    {
        // the atoms that have not changed in the new peptide
        ArrayList<Atom> oldAtoms = new ArrayList<>();

        // the atoms that are new in newPeptide
        ArrayList<Atom> newAtoms = new ArrayList<>();

        // populate lists
        for (Atom a : contents)
            {
                if ( oldMolecule.contents.contains(a) )
                    oldAtoms.add(a);
                else
                    newAtoms.add(a);
            }

        // compare distances between old atoms and new atoms
        for (Atom oldAtom : oldAtoms)
            {
                Vector3D oldPosition = oldAtom.position;
                for (Atom newAtom : newAtoms)
                    {
                        Vector3D newPosition = newAtom.position;
                        double distance = Vector3D.distance(oldPosition,newPosition);
                        if ( distance < Settings.MINIMUM_INTERATOMIC_DISTANCE &&
                             connectivity.getEdge(oldAtom, newAtom) == null )
                            return true;
                    }
            }
        return false;
    }

    /**
     * Checks if the atoms are too close in a molecule.  The minimum distance
     * is controlled by Settings.MINIMUM_DISTANCE.  Computations are performed
     * on the triangular distance matrix.
     * @return true if there is at least one atom that is too close to another atom
     */
    public boolean tooClose()
    {
        // compute the upper triangle of distance contacts
        for (int i=0; i < contents.size(); i++)
            {
                Atom atom1 = contents.get(i);
                Vector3D position1 = atom1.position;
                for (int j=i+1; j < contents.size(); j++)
                    {
                        Atom atom2 = contents.get(j);
                        Vector3D position2 = atom2.position;

                        // ignores distances between directly connected atoms
                        if ( Vector3D.distance(position1, position2) < Settings.MINIMUM_INTERATOMIC_DISTANCE &&
                             connectivity.getEdge(atom1,atom2) == null )
                            return true;
                    }
            }
        return false;
    }

    /**
     * Calculates the "distance" between two molecules by comparing their positions.
     * @param molecule1 the first molecule
     * @param molecule2 the second molecule
     * @return the RMS average distance between the atoms in molecule1 and molecule2 in angstroms
     */
    public static double calculateRMSD(Molecule molecule1, Molecule molecule2)
    {
        if ( molecule1.contents.size() != molecule2.contents.size() )
            throw new IllegalArgumentException("molecules are not the same size");
        
        double RMSD = 0.0;
        for (int i=0; i < molecule1.contents.size(); i++)
            {
                Vector3D fromVector = molecule1.contents.get(i).position;
                Vector3D toVector = molecule2.contents.get(i).position;
                RMSD += Math.pow(Vector3D.distance(fromVector,toVector), 2);
            }
        RMSD = RMSD * ( 1.0 / molecule1.contents.size() );
        RMSD = Math.sqrt(RMSD);
        return RMSD;
    }

    /**
     * Calculates the "distance" between two molecules by comparing their positions.
     * Molecules are superimposed first.
     * @param molecule1 the first molecule
     * @param molecule2 the second molecule
     * @param atomIndices the indices (0,1,...,n-1) of the atoms to be compared
     * @return the RMS average distance between the superimposed molecules, counting only the atoms in the lsit
     */
    public static double calculateRMSD(Molecule molecule1, Molecule molecule2, List<Integer> atomIndices)
    {
        Molecule m2 = superimpose(molecule1,molecule2,atomIndices);

        double RMSD = 0.0;
        for (Integer atomIndex : atomIndices)
            {
                Atom fromAtom = molecule1.contents.get(atomIndex);
                Atom toAtom = m2.contents.get(atomIndex);
                Vector3D fromVector = fromAtom.position;
                Vector3D toVector = toAtom.position;
                RMSD += Math.pow(Vector3D.distance(fromVector,toVector), 2);
            }
        RMSD = RMSD * ( 1.0 / molecule1.contents.size() );
        RMSD = Math.sqrt(RMSD);
        return RMSD;
    }

    /**
     * Superimposes two molecules based on a list of atom numbers.
     * Superposition is performed using the
     * <a href="http://en.wikipedia.org/wiki/Kabsch_algorithm">Kabsch algorithm</a>.
     * @param molecule1 one of the molecules
     * @param molecule2 another molecule (superposition transformation will be applied to this vector)
     * @param atomIndices the atom indices to use for the superposition
     * @return a new Molecule, which is molecule2 rotated to be maximally superimposed with molecule1
     */
    public static Molecule superimpose(Molecule molecule1, Molecule molecule2, List<Integer> atomIndices)
    {
        // check validity of molecule arguments
        if ( molecule1.contents.size() != molecule2.contents.size() )
            throw new IllegalArgumentException("molecule size mismatch");

        // check validity of atom numbers list
        // does not check the validity of the numbers themselves
        if ( atomIndices == null || atomIndices.size() < 3 || atomIndices.size() > molecule1.contents.size() )
            throw new IllegalArgumentException("invalid atom numbers list");

        // collect centroids
        Vector3D centroid1 = molecule1.getCentroid();
        Vector3D centroid2 = molecule2.getCentroid();

        // move molecules to origin
        Molecule m1 = molecule1.normalize();
        Molecule m2 = molecule2.normalize();

        // calculate superposition on centered molecules
        double[][] pre_matrixP = new double[atomIndices.size()][3];
        double[][] pre_matrixQ = new double[atomIndices.size()][3];
        for (int i=0; i < atomIndices.size(); i++)
            {
                int atomIndex = atomIndices.get(i);
                Atom fromAtom = m1.contents.get(atomIndex);
                if ( fromAtom == null )
                    throw new NullPointerException("error in atom index (" + atomIndex + ")");
                pre_matrixP[i][0] = fromAtom.position.getX();
                pre_matrixP[i][1] = fromAtom.position.getY();
                pre_matrixP[i][2] = fromAtom.position.getZ();

                Atom toAtom = m2.contents.get(atomIndex);
                pre_matrixQ[i][0] = toAtom.position.getX();
                pre_matrixQ[i][1] = toAtom.position.getY();
                pre_matrixQ[i][2] = toAtom.position.getZ();
            }

        Matrix matrixP = new Matrix(pre_matrixP);
        Matrix matrixQ = new Matrix(pre_matrixQ);
        Matrix matrixPtranspose = matrixP.transpose();

        // calculate covariance matrix and its singular value decomposition
        Matrix covarianceMatrix = matrixPtranspose.times(matrixQ);
        SingularValueDecomposition SVD = new SingularValueDecomposition(covarianceMatrix);
        Matrix matrixV = SVD.getU();
        Matrix matrixW = SVD.getV();

        // check sign of rotation matrix
        Double determinant = covarianceMatrix.det();
        double sign = 0.0;
        if (determinant > 0)
            sign = 1.0;
        else if (determinant == 0)
            sign = 1.0;
        else if (determinant < 0)
            sign = -1.0;

        double[][] pre_signedMatrix = {{1.0,0.0,0.0},{0.0,1.0,0.0},{0.0,0.0,sign}};
        Matrix signedMatrix = new Matrix(pre_signedMatrix);

        // compute optimal rotation matrix
        Matrix rotationMatrix = matrixV.times(signedMatrix);
        rotationMatrix = rotationMatrix.times(matrixW.transpose());

        // translate between packages by converting Vector3D objects to Matrix objects
        double[][] currentVectorPreMatrix = new double[3][1];
        double[][] newVectorPreMatrix = new double[3][1];
        
        // apply rotation to all atoms
        Map<Atom,Atom> atomMap = new HashMap<>();
        for (int i=0; i < m2.contents.size(); i++)
            {
                Atom fromAtom = m2.contents.get(i);
                Vector3D currentVector = fromAtom.position;
                currentVectorPreMatrix[0][0] = currentVector.getX();
                currentVectorPreMatrix[1][0] = currentVector.getY();
                currentVectorPreMatrix[2][0] = currentVector.getZ();
                Matrix currentVectorMatrix = new Matrix(currentVectorPreMatrix);

                Matrix newVectorMatrix = rotationMatrix.times(currentVectorMatrix);
                newVectorPreMatrix = newVectorMatrix.getArrayCopy();
                Vector3D newVector = new Vector3D(newVectorPreMatrix[0][0], newVectorPreMatrix[1][0], newVectorPreMatrix[2][0]);
                
                Atom newAtom = fromAtom.moveAtom(newVector);
                atomMap.put(fromAtom, newAtom);
            }

        // rotate molecule and undo translation
        Molecule rotatedMolecule = m2.moveAtoms(atomMap);
        return rotatedMolecule.shift(centroid2);
    }

    /**
     * Returns a deep copy of this molecule with a new name.
     * @param name the new name
     * @return the new molecule
     */
    public Molecule setName(String name)
    {
        Molecule copiedMolecule = moveAtoms(new HashMap<Atom,Atom>());
        return new Molecule(name, copiedMolecule.contents, copiedMolecule.connectivity);
    }
    
    @Override
    public int hashCode()
    {
        return Objects.hash(name, contents, connectivity);
    }

    /**
     * Checks for object equality by comparing all fields.
     * @param obj the object to be compared to
     * @return true if equal
     */
    @Override
    public boolean equals(Object obj)
    {
        if ( obj == null )
            return false;
        if ( obj == this )
            return true;
        if ( !(obj instanceof Molecule) )
            return false;

        Molecule anotherMolecule = (Molecule)obj;
        if ( this.name.equals(anotherMolecule.name) &&
             this.contents.equals(anotherMolecule.contents) &&
             this.connectivity.equals(anotherMolecule.connectivity) )
            return true;
        return false;
    }

    @Override
    public String toString()
    {
        String returnString = name + "\n\n";
        for (Atom a : contents)
            returnString = returnString + (contents.indexOf(a) + 1) + a.toString() + "\n";
        return returnString;
    }

    /**
     * Returns a GaussView readable string.
     */
    public String toGaussianString()
    {
        String geometry = "";
        for (Atom a : contents)
            geometry = geometry + a.toString() + "\n";
        return String.format("#\n\n%s\n\n0 1\n%s\n\n", name, geometry, connectivity);
    }

    /**
     * Returns the geometry string for use in Omnisol.
     */
    public String toOmnisolString()
    {
        String geometry = "";
        for (Atom a : contents)
            geometry = geometry + a.toString() + "\n";
        return geometry;
    }

    /**
     * Returns the string representation of a Molecule in XYZ files (to be passed to Tinker).
     * @param forcefield tells us which forcefield type to use
     * @return the xyz string
     */
    public String toXYZString(Forcefield forcefield) 
    {
	    //creates String to write
        String outputString = "";
        
        //write number of atoms and molecule name
        outputString = contents.size() + " " + name + "\n";

        //write atom list and connections
        int currentAtomNumber = 1;
        for (Atom currentAtom : contents)
            {
                Set<DefaultWeightedEdge> bonds = connectivity.edgesOf(currentAtom);
                Integer atomType = null;
                if ( forcefield == Forcefield.AMOEBA )
                    atomType = currentAtom.type1;
                else if ( forcefield == Forcefield.OPLS )
                    atomType = currentAtom.type2;
                else
                    throw new IllegalArgumentException("unknown forcefield type");
                outputString = outputString + String.format("%3d %2s %12.8f %12.8f %12.8f %6d", currentAtomNumber,
                                                            currentAtom.element.symbol,  currentAtom.position.getX(),
                                                            currentAtom.position.getY(), currentAtom.position.getZ(),
                                                            atomType);

                for (DefaultWeightedEdge e : bonds)
                    {
                        int edgeTarget = contents.indexOf(connectivity.getEdgeTarget(e)) + 1;
                        int edgeSource = contents.indexOf(connectivity.getEdgeSource(e)) + 1;
                        int edgeTargetToWrite = edgeTarget;
                        //change to edgeSource if edgeTarget is simply the current atom (edges have no particular ordering)
                        if (edgeTarget == currentAtomNumber)
                            edgeTargetToWrite = edgeSource;
                        outputString = outputString+String.format("%6d",edgeTargetToWrite);
                    }
                outputString = outputString + "\n";
                currentAtomNumber++;
            }
        return outputString;
    }
}
