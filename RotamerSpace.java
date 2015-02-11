import java.util.*;
import java.io.*;
import org.apache.commons.math3.geometry.euclidean.threed.*;
import org.jgrapht.*;
import org.jgrapht.graph.*;
import org.jgrapht.alg.*;
import com.google.common.collect.*;
import java.util.concurrent.*;
import java.util.concurrent.atomic.*;

/**
 * This represents the possible rotamers at each position of a peptide.
 */
public class RotamerSpace implements Immutable
{
    /** the peptide to pack */
    public final Peptide peptide;

    /**
     * the possible rotamers at each position of the sequence.
     * outer list is the position in the sequence.
     * inner list is a possible rotamer.
     * rotamers are the same amino acid/sidechain.
     */
    public final List<List<Rotamer>> rotamerSpace;

    /** pairs of incompatible rotamers */
    public final Set<RotamerPair> incompatiblePairs;

    /**
     * Creates a RotamerSpace just by copying fields.
     * @param peptide the peptide to pack the rotamers of
     * @param rotamerSpace the allowed rotamers at each position (an empty list means the position is fixed)
     * @param incompatiblePairs pairs of incompatible rotamers
     */
    public RotamerSpace(Peptide peptide, List<List<Rotamer>> rotamerSpace, Set<RotamerPair> incompatiblePairs)
    {
        this.peptide = peptide;
        this.rotamerSpace = rotamerSpace;
        this.incompatiblePairs = incompatiblePairs;
    }

    /**
     * Creates a RotamerSpace.  Assumes the inputPeptide is clash-free.  Throws an exception if no solutions are possible.
     * The input peptide will be mutated if there are any variable positions with only one possible rotamer.
     * @param inputPeptide the peptide to pack the rotamers of
     * @param whether HNs should be considered part of the rotamers
     */
    public static RotamerSpace create(Peptide inputPeptide, boolean includeHN)
    {
        // get rotamer space, which may or may not span different amino acids
        List<List<Rotamer>> rotamerSpace = getRotamerSpace(inputPeptide, includeHN);

        // if a position has only one possible rotamer, reconstitute the peptide with that rotamer
        Peptide peptide = inputPeptide;
        for (int i=0; i < rotamerSpace.size(); i++)
            {
                if ( rotamerSpace.get(i).size() == 1 )
                    {
                        Rotamer onlyRotamer = rotamerSpace.get(i).get(0);
                        peptide = Rotamer.reconstitute(peptide, onlyRotamer);
                    }
            }

        // determine which positions are variable (i.e., ones where we have a choice of rotamers)
        List<Integer> variablePositions = new ArrayList<>(peptide.sequence.size());
        for (int i=0; i < rotamerSpace.size(); i++)
            if ( rotamerSpace.get(i).size() > 1 )
                variablePositions.add(i);

        // get backbone atoms, which include the sidechains for fixed positions
        // outer index is residue index
        List<List<Atom>> backboneAtoms = getBackboneAtoms(peptide, variablePositions, includeHN);

        // prune rotamers that clash with the backbone
        rotamerSpace = pruneRotamerSpace(peptide, backboneAtoms, rotamerSpace);

        // check if any solutions are possible
        for (Integer i : variablePositions)
            {
                List<Rotamer> list = rotamerSpace.get(i);
                if ( list.size() == 0 )
                    throw new IllegalArgumentException("no rotamers left at position " + i + " after initial pruning");
            }

        // figure out which pairs are incompatible
        Set<RotamerPair> incompatiblePairs = Collections.newSetFromMap(new ConcurrentHashMap<RotamerPair,Boolean>());
        PairIterator iterator = new PairIterator(rotamerSpace, 1000);
        List<Future<Result>> futures = new ArrayList<>();
        while (iterator.hasNext())
            {
                List<RotamerPair> thisBatch = iterator.next();
                IncompatibleWorkUnit unit = new IncompatibleWorkUnit(rotamerSpace,thisBatch,incompatiblePairs);
                Future<Result> f = GeneralThreadService.submit(unit);
                futures.add(f);
            }

        GeneralThreadService.silentWaitForFutures(futures);
        
        // prune incompatible rotamers
        pruneIncompatibleRotamers(rotamerSpace, incompatiblePairs);
   
        // check if any solutions are possible
        checkRotamerSpace(rotamerSpace, peptide, incompatiblePairs);

        // return result
        printRotamerSizes(rotamerSpace);
        return new RotamerSpace(peptide, rotamerSpace, incompatiblePairs);
    }

    /**
     * Figures out which rotamers are possible for a given peptide.
     * @param inputPeptide the peptide to analyze
     * @param includeHN whether to consider backbone HNs part of sidechains
     * @return the rotamer space (empty inner lists mean no variation is desired)
     */
    public static List<List<Rotamer>> getRotamerSpace(Peptide inputPeptide, boolean includeHN)
    {
        throw new IllegalArgumentException("must override");
    }

    /**
     * Returns the atoms that are in the backbone.  Backbone atoms are defined as those that will never move during
     * rotamer packing.  If a position is marked as variable, its sidechain atoms will not be included in the result.
     * @param peptide the peptide to analyze
     * @param variablePositions indices in the peptide sequence where we want to vary the rotamers
     * @param includeHN if true, the backbone HNs are considered part of the sidechain
     * @return the atoms in the backbone, indexed by residue
     */
    public static List<List<Atom>> getBackboneAtoms(Peptide peptide, List<Integer> variablePositions, boolean includeHN)
    {
        List<List<Atom>> backboneAtoms = new ArrayList<>(peptide.sequence.size());
        for (int i=0; i < peptide.sequence.size(); i++)
            {
                Residue residue = peptide.sequence.get(i);
                List<Atom> atoms = new ArrayList<>(residue.atoms);
                
                // if the position is variable, do not add the sidechain atoms
                if ( variablePositions.contains(i) )
                    atoms.removeAll( getSidechainAtoms(peptide, residue, includeHN) );
                
                backboneAtoms.add(ImmutableList.copyOf(atoms));
           }
        return ImmutableList.copyOf(backboneAtoms);
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
                sidechainAtoms = RotamerFactory.getHalfGraph(connectivity, atomA, atomB);
            }
        else
            {
                // just add the sidechain atoms
                sidechainAtoms = peptide.getHalfGraph(atomA, atomB);
            }
        return sidechainAtoms;
    }

    /**
     * Takes the current rotamer space and removes those that clash with the backbone.
     * Resulting lists are immutable.
     * @param startingPeptide the structure the backbone atoms are from
     * @param backboneAtoms the backbone atoms
     * @param unprunedRotamerSpace the various rotamers possible at each position
     * @return a list of lists containing the pruned rotamers
     */
    public static List<List<Rotamer>> pruneRotamerSpace(Peptide startingPeptide, List<List<Atom>> backboneAtoms, List<List<Rotamer>> unprunedRotamerSpace)
    {
        List<List<Rotamer>> returnList = new ArrayList<>(unprunedRotamerSpace.size());
        for (int i=0; i < startingPeptide.sequence.size(); i++)
            {
                List<Rotamer> unprunedList = unprunedRotamerSpace.get(i);
                List<Rotamer> prunedList = pruneRotamers(startingPeptide, backboneAtoms, unprunedList);
                returnList.add(prunedList);
            }
        return ImmutableList.copyOf(returnList);
    }

    /**
     * Prunes any rotamers from the given list that clash with the backbone atoms.
     * Resulting list is immutable.
     * @param backboneAtoms the backbone atoms
     * @param thisRotamerSpace the various rotamers possible at this particular position in the sequence
     * @return the pruned list
     */
    public static List<Rotamer> pruneRotamers(Peptide startingPeptide, List<List<Atom>> backboneAtoms, List<Rotamer> thisRotamerSpace)
    {
        List<Rotamer> returnList = new LinkedList<>();
        for (Rotamer rotamer : thisRotamerSpace)
            {
                boolean clashes = checkRotamerBackboneClash(startingPeptide, backboneAtoms, rotamer);
                if ( !clashes )
                    returnList.add(rotamer);
            }
        return ImmutableList.copyOf(returnList);
    }

    /**
     * Checks if the specified rotamer clashes with the backbone.  Assumes the rotamer will not clash
     * with itself.  Note that the rotamer might include atoms that aren't in the peptide.
     * @param startingPeptide the structure the backbone atoms are from
     * @param backboneAtoms the backbone atoms
     * @param rotamer the rotamer to check
     * @return true if the rotamer clashes
     */
    public static boolean checkRotamerBackboneClash(Peptide startingPeptide, List<List<Atom>> backboneAtoms, Rotamer rotamer)
    {
        // if this is glycine, automatically let it go by
        if ( rotamer.description.indexOf("glycine") > -1 )
            return false;

        // if this is a trans-proline, check the improper torsion angle at the proline
        // and reject it if it's too non-planar
        Residue residue = startingPeptide.sequence.get(rotamer.sequenceIndex);
        if ( residue.aminoAcid.isProline() && Math.abs(residue.omega.getDihedralAngle() - 180.0) < 30.0 )
            {
                Atom Cdelta_i = null;
                for (Atom a : rotamer.atoms)
                    {
                        if ( a.type1 == 59 )
                            {
                                Cdelta_i = a;
                                break;
                            }
                    }
                Atom Calpha_i = residue.phi.atom3;
                Atom N_i = residue.phi.atom2;
                Atom C_i_minus1 = residue.phi.atom1;

                ProtoTorsion improperTorsion = new ProtoTorsion(C_i_minus1, N_i, Calpha_i, Cdelta_i);
                double angle = improperTorsion.getDihedralAngle();
                //System.out.printf("improper: %.1f\n", angle);
                
                // don't allow really non-planar nitrogen atoms
                if ( Math.abs(angle-180.0) > 30.0 )
                    {
                        //System.out.println("rejected");
                        return true;
                    }
            }

        // if this is a proline, remove the HN of the current residue from consideration
        // otherwise, there will be a hydrogen right on top of the delta atoms of the proline sidechain
        // and we'll artificially get a clash
        Atom prolineHN = null;
        if ( residue.aminoAcid.isProline() )
            prolineHN = residue.HN;

        // check for clashes
        for (int i=0; i < backboneAtoms.size(); i++)
            {
                // don't check clashes with the rotamer's own backbone atoms
                if ( i == rotamer.sequenceIndex )
                    continue;
                List<Atom> backboneAtomsInThisResidue = backboneAtoms.get(i);
                for (Atom backboneAtom : backboneAtomsInThisResidue)
                    {
                        if ( backboneAtom == prolineHN )
                            {
                                //System.out.println("skipped backbone HN");
                                continue;
                            }
                        
                        for (Atom rotamerAtom : rotamer.atoms)
                            {
                                double distance = Molecule.getDistance(backboneAtom, rotamerAtom);
                                if ( distance < Settings.MINIMUM_INTERATOMIC_DISTANCE )
                                    return true;
                            }
                    }
            }
        return false;
    }

    /**
     * Counts the total number of rotamers in a rotamer space.
     * @param rotamerSpace the rotamer space to count up
     * @return the number of rotamers
     */
    public static int countRotamers(List<List<Rotamer>> rotamerSpace)
    {
        int count = 0;
        for (List<Rotamer> list : rotamerSpace)
            for (Rotamer r : list)
                count++;
        return count;
    }

    /**
     * Prints out how many rotamers there are at each position on one line
     * @param rotamerSpace the rotamer space to print out
     */
    public static void printRotamerSizes(List<List<Rotamer>> rotamerSpace)
    {
        for (List<Rotamer> list : rotamerSpace)
            System.out.print(list.size() + " ");
        System.out.println(" / " + countRotamers(rotamerSpace) + " total");
    }
       
    /**
     * Removes rotamers that cannot possibly be part of the global minimum energy conformation (GMEC).
     * Mutates the parameters.  By definition, rotamers that are incompatible with rotamers in the GMEC
     * cannot themselves be part of the GMEC.  If a rotamer is the only possibility at a position, then it is
     * by definition a part of the GMEC.  This method will also prune the incompatible rotamer pairs.
     * @param rotamerSpace the rotamer space (will be modified)
     * @param incompatiblePairs the incompatible pairs (will be modified)
     */
    public static void pruneIncompatibleRotamers(List<List<Rotamer>> rotamerSpace, Set<RotamerPair> incompatiblePairs)
    {
        // first, find the GMEC rotamers
        HashSet<Rotamer> GMECrotamers = new HashSet<>();
        for (List<Rotamer> list : rotamerSpace)
            {
                if ( list.size() == 1 )
                    GMECrotamers.add(list.get(0));
            }

        // second, locate incompatible pairs containing GMEC rotamers and load their
        // partners into a set of things to be removed
        HashSet<Rotamer> rotamersToBeRemoved = new HashSet<>();
        HashSet<RotamerPair> rotamerPairsToBeRemoved = new HashSet<>();
        for (RotamerPair pair : incompatiblePairs)
            {
                Rotamer rotamer1 = pair.rotamer1;
                Rotamer rotamer2 = pair.rotamer2;
                boolean contains1 = GMECrotamers.contains(rotamer1);
                boolean contains2 = GMECrotamers.contains(rotamer2);
                if ( contains1 && contains2 )
                    throw new IllegalArgumentException("two rotamers in the GMEC cannot be incompatible!");
                else if ( contains1 )
                    rotamersToBeRemoved.add(rotamer2);
                else if ( contains2 )
                    rotamersToBeRemoved.add(rotamer1);

                if ( contains1 || contains2 )
                    rotamerPairsToBeRemoved.add(pair);
            }

        // remove the inadmissible rotamers
        for (List<Rotamer> list : rotamerSpace)
            list.removeAll(rotamersToBeRemoved);

        // delete the unnecessary incompatible pairs
        incompatiblePairs.removeAll(rotamerPairsToBeRemoved);
        //if ( rotamersToBeRemoved.size() > 0 || rotamerPairsToBeRemoved.size() > 0 )
        //    System.out.printf("%d unnecessary rotamers and %d incompatible pairs have been removed.\n", rotamersToBeRemoved.size(), rotamerPairsToBeRemoved.size());
    }

    /**
     * Ensures that it is at least superficially possible to find a solution.
     * Throws an exception if the rotamer space is invalid.
     * @param rotamerSpace the rotamer space
     * @param peptide the parent peptide
     * @param incompatiblePairs pairs of incompatible rotamers - optional
     */
    public static void checkRotamerSpace(List<List<Rotamer>> rotamerSpace, Peptide peptide, Set<RotamerPair> incompatiblePairs)
    {
        // if we prune all the rotamers from a non-hairpin position, abort
        int sequenceLength = peptide.sequence.size();
        int forbiddenIndex = (sequenceLength/2) - 1;
        for (int i=0; i < rotamerSpace.size(); i++)
            {
                if ( i == forbiddenIndex || i == forbiddenIndex + 1 || rotamerSpace.get(i).size() > 0 )
                    continue;
                throw new IllegalArgumentException(String.format("can't continue because there are no valid rotamers at a non-hairpin position (%d)", i));
            }

        // ensure that every position can make at least one compatible pair with the rotamers at all other positions
        // doesn't guarantee there's a solution
        if ( incompatiblePairs == null)
            throw new NullPointerException("null incompatible pairs");
        for (int i=0; i < rotamerSpace.size()-1; i++)
            {
                List<Rotamer> list1 = rotamerSpace.get(i);
                // skip over empty positions
                if ( list1.size() == 0 )
                    continue;

                for (int j=i+1; j < rotamerSpace.size(); j++)
                    {
                        List<Rotamer> list2 = rotamerSpace.get(j);
                        if ( list2.size() == 0 )
                            continue;

                        // try to find at least one commpatible match
                        boolean positionOK = false;
                        rotamers:
                        for (Rotamer rotamer1 : list1)
                            {
                                for (Rotamer rotamer2 : list2)
                                    {
                                        RotamerPair pair = new RotamerPair(rotamer1, rotamer2);
                                        if ( !incompatiblePairs.contains(pair) )
                                            {
                                                positionOK = true;
                                                break rotamers;
                                            }
                                    }
                            }
                        if ( !positionOK )
                            throw new IllegalArgumentException("There are no valid solutions for position " + i + ".");
                    }
            }

        //System.out.println("Rotamer space passes preliminary checks.");
    }

    /**
     * Checks a batch of rotamer pairs and figures out which ones are incompatible.
     */
    public static class IncompatibleWorkUnit implements WorkUnit
    {
        /** the rotamer space to search */
        private final List<List<Rotamer>> rotamerSpace;

        /** the pairs of rotamers to compare */
        private final List<RotamerPair> work;

        /** where to put the results */
        private final Set<RotamerPair> incompatiblePairs;

        public static AtomicInteger counter = new AtomicInteger();
        public final int count;

        /** constructor */
        public IncompatibleWorkUnit(List<List<Rotamer>> rotamerSpace, List<RotamerPair> work, Set<RotamerPair> incompatiblePairs)
        {
            this.rotamerSpace = rotamerSpace;
            this.work = work;
            this.incompatiblePairs = incompatiblePairs;
            this.count = counter.getAndIncrement();
        }

        /** analyzes the pairs of rotamers */
        public Result call()
        {
            // analyze the rotamer pairs
            HashSet<RotamerPair> thisIncompatible = new HashSet<>();
            for (RotamerPair pair : work)
                {
                    if ( isIncompatible(pair) ) 
                        thisIncompatible.add(pair);
                }
            
            // incorporate results
            incompatiblePairs.addAll(thisIncompatible);
            return null;
        }

        /**
         * Determines whether two rotamers are incompatible.
         * @param pair the pair of rotamers to check
         * @return true if the rotamers are incompatible
         */
        public static boolean isIncompatible(RotamerPair pair)
        {
            for (Atom a1 : pair.rotamer1.atoms)
                {
                    for (Atom a2 : pair.rotamer2.atoms)
                        {
                            double distance = Vector3D.distance(a1.position, a2.position);
                            if ( distance < Settings.MINIMUM_INTERATOMIC_DISTANCE )
                                return true;
                        }
                }
            return false;
        }

        @Override
        public String toString()
        {
            return "IncompatibleWorkUnit with " + work.size() + " pairs of rotamers";
        }
    }

    /**
     * Given a rotamer space, create all the possible pairs of rotamers without
     * any duplicates or intra-residue pairs.
     */
    public static class PairIterator
    {
        /** the rotamer space to traverse */
        private final List<List<Rotamer>> rotamerSpace;

        /** whether the iterator has anything to give next */
        private boolean hasNext;

        /** how many rotamer pairs to return at a time */
        private final int batchSize;

        /**
         * the current pair we are producing
         */
        private IndexPair currentPosition;

        /** constructor */
        public PairIterator(List<List<Rotamer>> rotamerSpace, int batchSize)
        {
            this.rotamerSpace = rotamerSpace;
            this.batchSize = batchSize;

            // find the first residue that has a non-zero number of rotamers
            int index1 = -1;
            int index3 = -1;
            for (List<Rotamer> list : rotamerSpace)
                {
                    if ( index1 == -1 && list.size() > 0 )
                        index1 = rotamerSpace.indexOf(list);
                    else if ( index3 == -1 && list.size() > 0 )
                        index3 = rotamerSpace.indexOf(list);
                }
            if ( index1 == -1 || index3 == -1 )
                throw new IllegalArgumentException("nothing to start the iterator with");
            hasNext = true;
            currentPosition = new IndexPair(index1,0,index3,0);
        }

        /**
         * Returns the rotamer pair corresponding to these indices.
         */
        private RotamerPair getRotamerPair(int index1, int index2, int index3, int index4)
        {
            Rotamer rotamer1 = rotamerSpace.get(index1).get(index2);
            Rotamer rotamer2 = rotamerSpace.get(index3).get(index4);
            return new RotamerPair(rotamer1, rotamer2);
        }

        /**
         * Returns the next list we can use as the source for a rotamer pair.
         * @param currentFromList the current source
         * @return the next list
         */
        private List<Rotamer> nextFromList(List<Rotamer> currentFromList)
        {
            List<Rotamer> nextFromList = null;
            for (int i=rotamerSpace.indexOf(currentFromList)+1; i < rotamerSpace.size() - 1; i++)
                {
                    List<Rotamer> list = rotamerSpace.get(i);
                    if ( list.size() > 0 )
                        {
                            nextFromList = list;
                            break;
                        }
                }
            return nextFromList;
        }

        /**
         * Returns the next list we can use as the target for a rotamer pair.
         * @param currentToList the current target
         * @return the next list
         */
        private List<Rotamer> nextToList(List<Rotamer> currentToList)
        {
            List<Rotamer> nextToList = null;
            for (int i=rotamerSpace.indexOf(currentToList)+1; i < rotamerSpace.size(); i++)
                {
                    List<Rotamer> list = rotamerSpace.get(i);
                    if ( list.size() > 0 )
                        {
                            nextToList = list;
                            break;
                        }
                }
            return nextToList;
        }

        /**
         * Whether we are finished with this iterator.
         * @return false if we are done
         */
        public boolean hasNext()
        {
            return hasNext;
        }

        /**
         * Return the next batch of rotamer pairs.  The currentPosition will be set at
         * the first thing of the next map to return on the next call of next().
         * @return the next batch of rotamer pairs
         */
        public List<RotamerPair> next()
        {
            // if we are finished, throw an exception
            if ( !hasNext )
                throw new IllegalArgumentException("iterator is finished");
            // otherwise, make the next batch
            
            // the current position
            int index1 = currentPosition.outerIndex1;
            int index2 = currentPosition.innerIndex1;
            int index3 = currentPosition.outerIndex2;
            int index4 = currentPosition.innerIndex2;

            // create the next batch
            List<Rotamer> fromList = rotamerSpace.get(index1);
            Rotamer fromRotamer = fromList.get(index2);
            List<Rotamer> toList = rotamerSpace.get(index3);
            Rotamer toRotamer = toList.get(index4);

            // iterate through the current batch
            ArrayList<RotamerPair> nextBatch = new ArrayList<>(batchSize);
            int size = 0;
            while ( size < batchSize )
                {
                    // add the current pair
                    nextBatch.add(new RotamerPair(fromRotamer, toRotamer));
                    size++;

                    // increment the current position
                    index4++;

                    // if we have run out for the current to list
                    if ( index4 > toList.size() - 1 )
                        {
                            // try to move onto the next to list
                            index4 = 0;
                            toList = nextToList(toList);

                            // if we have run out of to lists
                            if ( toList == null )
                                {
                                    // try to move onto the next entry in the from list
                                    index2++;

                                    // if we have run out for the current from list
                                    if ( index2 > fromList.size() - 1 )
                                        {
                                            // we have run out of rotamers in this from list, so
                                            // attempt to move to the next from list
                                            index2 = 0;
                                            fromList = nextFromList(fromList);

                                            if ( fromList != null )
                                                {
                                                    // get the first rotamer
                                                    fromRotamer = fromList.get(index2);
                                                    toList = nextToList(fromList);
                                                    if ( toList == null )
                                                        {
                                                            hasNext = false;
                                                            return nextBatch;
                                                        }
                                                }
                                            else
                                                {
                                                    // we have run out of from lists, so we are finished iterating
                                                    hasNext = false;
                                                    return nextBatch;
                                                }
                                        }
                                    else
                                        {
                                            // entries are still available in the from list, so
                                            // get the new from rotamer
                                            fromRotamer = fromList.get(index2);
                                            toList = nextToList(fromList);
                                            if ( toList == null )
                                                {
                                                    hasNext = false;
                                                    return nextBatch;
                                                }
                                        }
                                }
                        }

                    // set the next pair
                    toRotamer = toList.get(index4);
                }

            // set the current position
            index1 = rotamerSpace.indexOf(fromList);
            index3 = rotamerSpace.indexOf(toList);
            currentPosition = new IndexPair(index1,index2,index3,index4);

            // return the current next
            return nextBatch;
        }
    }
}