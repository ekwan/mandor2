import java.util.*;
import java.io.*;
import org.apache.commons.math3.geometry.euclidean.threed.*;
import org.jgrapht.*;
import org.jgrapht.graph.*;
import org.jgrapht.alg.*;
import com.google.common.collect.*;
import java.util.concurrent.*;
import java.util.concurrent.atomic.*;
import com.google.common.util.concurrent.AtomicDouble;

/**
 * This class represents a memory-limited, parallel A* search for the lowest energy peptide given a
 * rotamer space.  The algorithm is described in doi:10.1093/bioinformatics/btu264.
 */
public class RotamerIterator
{
    // static fields

    /** the number of worker threads to use */
    public static final int NUMBER_OF_QUEUES = 512;

    /** maximum number of seconds to iterate for */
    public static final double MAX_TIME = 180.0;

    /** maximum number of poses to return */
    public final int MAX_POSES;

    // member fields
    
    /** the pruned rotamer space to iterate over */
    public final List<List<Rotamer>> rotamerSpace;

    /** the self-energies */
    public final Map<Rotamer,Double> rotamerSelfEnergies;

    /** the interaction energies */
    public final Map<RotamerPair,Double> rotamerInteractionEnergies;

    /** energy minima for heuristic function calculation.  */
    public final Map<Rotamer,double[]> energyMinima = new ConcurrentHashMap<>();

    /** the root of the tree */
    public final Node root;

    /** the new nodes found on every iteration */
    public final Set<Node> newNodes;

    /** the sequence positions at which we expect to have rotamers in the goal node */
    public final ImmutableSet<Integer> targetIndices;

    /** where the answers are kept */
    public final List<Node> solutions;

    /** whether to do things in parallel */
    public final boolean parallelize;

    /** parallelization is turned on by default */
    public RotamerIterator(List<List<Rotamer>> rotamerSpace, Map<Rotamer,Double> rotamerSelfEnergies,
                           Map<RotamerPair,Double> rotamerInteractionEnergies, int MAX_POSES)
    {
        this(rotamerSpace, rotamerSelfEnergies, rotamerInteractionEnergies, MAX_POSES, true);
    }

    /**
     * Public constructor.
     * @param rotamerSpace the pruned rotamer space to iterate over
     * @param rotamerSelfEnergies the energies of the individual rotamers, including rotamer-backbone interactions
     * @param rotamerInteractionEnergies the interaction energies (it is assumed that the interaction energies will be null if incompatible)
     * @param MAX_POSES the maximum number of poses to return
     * @param parallelize whether to parallelize the calculations
     */
    public RotamerIterator(List<List<Rotamer>> rotamerSpace, Map<Rotamer,Double> rotamerSelfEnergies,
                           Map<RotamerPair,Double> rotamerInteractionEnergies, int MAX_POSES, boolean parallelize)
    {
        // copy fields
        this.rotamerSpace = rotamerSpace;
        this.rotamerSelfEnergies = rotamerSelfEnergies;
        this.rotamerInteractionEnergies = rotamerInteractionEnergies;
        this.MAX_POSES = MAX_POSES;
        this.parallelize = parallelize;

        // make a common set of newly explored nodes to be distributed
        newNodes = Collections.newSetFromMap(new ConcurrentHashMap<Node,Boolean>());

        // precompute some needed energies
        //System.out.println("Precomputing energy minima...");
        //RotamerSpace.printRotamerSizes("debug", rotamerSpace);
        computeEnergyMinima();
        //System.out.println("\nPrecomputation is complete.");

        // create the parent node and set the target indices
        List<Rotamer> rotamers = new ArrayList<>();
        Set<Integer> tempSet = new HashSet<>(); // for the target indices 
        for (int i=0; i < rotamerSpace.size(); i++)
            {
                List<Rotamer> list = rotamerSpace.get(i);
                if ( list.size() == 1 )
                    rotamers.add(list.get(0));
                else if ( list.size() > 1 ) 
                    tempSet.add(i);
            }
        root = new Node(rotamers);
        targetIndices = ImmutableSet.copyOf(tempSet);

        // keep the answers in here
        solutions = Collections.synchronizedList(new ArrayList<Node>());
    }

    /**
     * Produces the energy minima.  This is defined as min on u of E(i_r, k_u), where k != i.
     */
    public void computeEnergyMinima()
    {
        List<Future<Result>> futures = new ArrayList<>();
        for (int i = 0; i < rotamerSpace.size(); i++) {
            if ( rotamerSpace.get(i).size() == 1 )
                continue;
            for (Rotamer r : rotamerSpace.get(i)) {
                if ( parallelize )
                    futures.add(GeneralThreadService.submit(new MinJob(r)));
                else
                    new MinJob(r).call();
            }
        }
        GeneralThreadService.silentWaitForFutures(futures);
    }

    /**
     * Produces the lowest energy poses, if any.  Executes in parallel.  Nothing else should be running in parallel.
     */
    public List<Node> iterate()
    {
        // create the queues
        List<PriorityQueue<Node>> queues = new ArrayList<>(NUMBER_OF_QUEUES);
        int numberOfQueues = NUMBER_OF_QUEUES;
        if ( !parallelize )
            numberOfQueues = 1;
        for (int i=0; i < numberOfQueues; i++)
            queues.add( new PriorityQueue<Node>() );

        // load the root node into one queue
        queues.get(0).add(root);

        // while there are still threads with non-empty queues, and we haven't found the maximum number of
        // solutions or exceeded the maximum amount of time
        //System.out.println("Starting graph search.");
        List<List<Node>> batches = new ArrayList<>();
        for (int i=0; i < numberOfQueues; i++)
            batches.add(new ArrayList<Node>());
        Date start = new Date();
        
        // iterate
        int iteration = 0;
        while (true)
        {
            iteration++;

            // create the latch
            CountDownLatch latch = new CountDownLatch(numberOfQueues);

            // submit all jobs
            for (int i=0; i < numberOfQueues; i++)
                {
                    PriorityQueue<Node> queue = queues.get(i);
                    SearchUnit job = new SearchUnit(queue, latch, batches.get(i));
                    if ( parallelize )
                        GeneralThreadService.submit(job);
                    else
                        job.call();
                }

            // await the latch
            try
                {
                    latch.await();
                }
            catch (Exception e)
                {
                    e.printStackTrace();
                    break;
                }
            //System.out.printf("All nodes expanded.  %d solutions have been found.\n", solutions.size());

            // check if we should stop
            // if we should stop, signal all threads to stop
            Date now = new Date();
            double elapsedTime = ( now.getTime() - start.getTime() ) / 1000.0;
            if ( solutions.size() >= MAX_POSES )
                {
                    //System.out.printf("\nTermination condition reached: %d solutions found.\n", solutions.size());
                    break;
                }
            else if ( elapsedTime > MAX_TIME )
                {
                    //System.out.println("\nTermination condition reached: maximum time exceeded.");
                    break;
                }

            // parcel out the newly discovered nodes to all the queues
            int newNodesSize = newNodes.size();
            batches = redistribute(newNodes);
            int totalNodes = 0;
            int minUndetermined = 20;
            double minCost = 0.0;
            double maxCost = 0.0;
            boolean first = true;
            for (PriorityQueue<Node> queue : queues)
                {
                    totalNodes += queue.size();
                    /*if ( queue.size() > 0 )
                        {
                            minUndetermined = Math.min(minUndetermined, queue.peek().undetermined.size());
                            double cost = queue.peek().totalCost;
                            if ( first || cost < minCost )
                                minCost = cost;
                            if ( first || cost > maxCost )
                                maxCost = cost;
                            first = false;
                        }*/
                }
            /*double throughput = totalNodes / elapsedTime;
            Date now2 = new Date();
            long overhead = now2.getTime() - now.getTime();
            //System.out.printf("iter %-5d   %3d solns   %5d new nodes    %2d min_undet    %8.2e tot_nodes   %8.1f minCost   %8.1f maxCost  %3d overhead   %7.1f elapsed   %7.1f nodes/s\r",
            //                  iteration, solutions.size(), newNodesSize, minUndetermined, totalNodes/1.0, minCost, maxCost, overhead, elapsedTime, throughput);
            */

            if ( totalNodes == 0 && newNodesSize == 0 )
                {
                    //System.out.println("\nNo nodes left.");
                    break;
                }
        }

        // sort and return the solutions
        Collections.sort(solutions);
        return solutions;
    }

    public class MinJob implements WorkUnit
    {
        /** the rotamer to minimize.  */
        public final Rotamer r;

        /** Constructor */
        public MinJob(Rotamer r)
        {
            this.r = r;
        }

        /** add the energy to energyMinima. */
        public Result call()
        {
            double[] output = new double[rotamerSpace.size()];
            for (int i = 0; i < output.length; i++) {
                if (i==r.sequenceIndex) continue;
                boolean first = true;
                for (Rotamer s : rotamerSpace.get(i)) {
                    Double e = rotamerInteractionEnergies.get(new RotamerPair(r,s));
                    if (e==null) continue; 
                    if (first) {
                        output[i] = e;
                        first = false;
                    } else if (e < output[i]) {
                        output[i] = e;
                    }
                }
            }
            // add the energies to energyMinima
            energyMinima.put(r,output);
            return null;
        }
    } // end of class MinJob

    /**
     * Redistributes the given nodes to the specified threads.  Could be used to redistributed newly discovered nodes or
     * re-parcel out nodes when queues get full.
     * @param newNodes the nodes to distribute (will be cleared after this operation)
     */
    public static List<List<Node>> redistribute(Set<Node> newNodes)
    {
        // don't bother if there's no work
//        if ( newNodes.size() == 0 )
//            return;
/*        System.out.println("Redistributing; newNodes.size() > 0.");

        // create a list of nodes where nodes with the same parents are contiguous
        List<Node> sorted = new LinkedList<>();
        while (newNodes.size() > 0)
            {
                System.out.println("newNodes.size() : " + newNodes.size());
                Node newNode = null;
                for (Node n : newNodes)
                    {
                        newNode = n;
                        break;
                    }
//                if (newNode!=null) {
                    boolean added = false;
                    for (int i=0; i < sorted.size()-1; i++)
                        {
                            Node currentNode = sorted.get(i);
                            if ( currentNode.sameParent(newNode) )
                                {
                                    sorted.add(i+1, newNode);
                                    added = true;
                                    System.out.println("newNode inserted at " + (i+1));
                                    break;
                                }
                        }
                    if ( !added ) {
                        sorted.add(newNode);
                        System.out.println("Adding Node to end of queue.");
                    }
                    newNodes.remove(newNode);
//                }
            }
*/

        
        // distribute the nodes
        List<List<Node>> batches = new ArrayList<>();
        for ( int i=0; i < NUMBER_OF_QUEUES; i++ )
            batches.add(new ArrayList<Node>());

        if ( newNodes.size() == 0 )
            return batches;

        int index = 0;
        for (Node node : newNodes)
            {
                if ( index > NUMBER_OF_QUEUES - 1 )
                    index = 0;
                batches.get(index).add(node);
                index++;
            }
        newNodes.clear();
        return batches;
    } // end of redistribute

    /** Represents a node on the rotamer tree. */
    public class Node implements Comparable<Node>
    {
        public final List<Rotamer> rotamers;

        /** the indices at which the rotamer is not yet determined.  */
        public final Set<Integer> undetermined;

        /** the value of g(x), the cost to get to this node */
        public final double actualCost;

        /** the value of f(x), which is g(x) + h(x), where h(x) is the heuristic distance to the goal */
        public final double totalCost;

        /** to instantiate multiple related nodes at the same time */
        private Node(List<Rotamer> rotamers, Set<Integer> undetermined, double actualCost, double totalCost)
        {
            this.rotamers = rotamers;
            this.undetermined = undetermined;
            this.actualCost = actualCost;
            this.totalCost = totalCost;
        }

        /** to initialize the parent node */
        public Node(List<Rotamer> rotamers)
        {
            this.rotamers = rotamers;
            this.undetermined = new TreeSet<>();
            for (int i = 0; i < rotamerSpace.size(); i++)
                if (rotamerSpace.get(i).size() > 0) undetermined.add(i);
            for (Rotamer r : rotamers) undetermined.remove(r.sequenceIndex);
        
            // calculate costs
            double tempActualCost = 0.0;
            for (Rotamer r : rotamers) tempActualCost += rotamerSelfEnergies.get(r);
            for (int i=0; i < rotamers.size()-1; i++)
                {
                    Rotamer rotamer1 = rotamers.get(i);
                    for (int j=i+1; j < rotamers.size(); j++)
                        {
                            Rotamer rotamer2 = rotamers.get(j);
                            RotamerPair pair = new RotamerPair(rotamer1, rotamer2);
                            tempActualCost += rotamerInteractionEnergies.get(pair);
                        }
                }
            actualCost = tempActualCost;
            totalCost = actualCost + heuristicCost();
        }

        /** to intialize descendents */
        public Node(Rotamer rotamer, Node parent)
        {
            rotamers = new ArrayList<>(parent.rotamers);
            rotamers.add(rotamer);
            undetermined = new TreeSet<>(parent.undetermined);
            boolean success = undetermined.remove(rotamer.sequenceIndex);
            if ( !success )
                throw new IllegalArgumentException("error removing");
            actualCost = actualCost(rotamer, parent);
            totalCost = actualCost + heuristicCost();
        }

        /**
         * Calculates the actual cost f(x) of this node, given the actual cost already
         * calculated in the parent node.  Only calculates the new terms to add to the actual cost
         * and reuses the old terms from the parent.
         * @param rotamer the rotamer we have just added to this node
         * @param parent the node that led to this one
         * @return the cost of this node
         */
        public double actualCost(Rotamer rotamer, Node parent)
        {
            double output = parent.actualCost;
            output += rotamerSelfEnergies.get(rotamer);
            for (Rotamer r : parent.rotamers)
                output += rotamerInteractionEnergies.get(new RotamerPair(rotamer,r));
            return output;
        }

        /**
         * Calculates the heuristic cost g(x) of this node, given a set of rotamers.  The
         * values of the minimum operator are gathered from a cache.
         * @return the heuristic cost
         */
        public double heuristicCost()
        {
            double output = 0.0;
            for (Integer i : undetermined) {
                double min = 0.0;
                boolean first = true;
                r:
                for (Rotamer r : rotamerSpace.get(i)) {
                    // self-energy
                    double energy = rotamerSelfEnergies.get(r);
                    // interaction energies with determined rotamers
                    for (Rotamer s : rotamers) {
                        Double e = rotamerInteractionEnergies.get(new RotamerPair(r,s));
                        if (e==null) continue r;
                        energy += e;
                    }
                    // interaction energies with undetermined rotamers
                    double[] minima = energyMinima.get(r);
                    for (Integer j : undetermined) if (j!=i) energy += minima[j];
                    // replace the minimum if necessary
                    if (first) {
                        min = energy;
                        first = false;
                    } else if (energy < min) {
                        min = energy;
                    }
                }
                output += min;
            }
            return output;
        }

        /**
         * Checks if there is a rotamer at all non-hairpin positions.
         * @return true if this is a goal node
         */
        public boolean isTarget()
        {
            return ( undetermined.size() == 0 );
        }

        /**
         * Checks if this has the same parent as another Node.  
         * @param other The other Node.  
         * @return true if these have the same parent.  
         */
        public boolean sameParent(Node other)
        {
            if (other==null) return false;
            if (this.rotamers.size()!=other.rotamers.size()) return false;
            for (int i = 0; i < this.rotamers.size()-1; i++)
                if (!this.rotamers.get(i).equals(other.rotamers.get(i))) return false;
            return true;
        }

        /**
         * Produce a list of all child Nodes of this Node.  Avoids duplicate work by sharing
         * some data and the undetermined field. 
         * @return A List of all children of this Node.  
         */
        public List<Node> expand()
        {
            // this will store the output
            List<Node> output = new ArrayList<>();

            // locate the next position to populate rotamers at
            int next = -1;
            for (Integer i : undetermined) {
                next = i;
                break;
            }

            // trying to expand a completed node is bad news
            if (next==-1)
                {
                    System.out.println("Warning: expanding an empty node.");
                    return new ArrayList<Node>();
                }

            // find rotamers that are compatible with those already present
            // store them in a list (newRotamers) and also construct the new lists of rotamers that
            // will go into the children (stored in newRotamerLists)
            List<Rotamer> newRotamers = new ArrayList<>(rotamerSpace.get(next).size());
            List<List<Rotamer>> newRotamerLists = new ArrayList<>(rotamerSpace.get(next).size());
            s:
            for (Rotamer s : rotamerSpace.get(next)) {
                for (Rotamer t : rotamers) {
                    Double e = rotamerInteractionEnergies.get(new RotamerPair(s,t));
                    if (e==null) continue s;
                }
                newRotamers.add(s);
                ArrayList<Rotamer> newList = new ArrayList<>(rotamers);
                newList.add(s);
                newRotamerLists.add(newList);
            }

            // make a new undetermined set that will be common to all children
            TreeSet<Integer> newUndetermined = new TreeSet<>(undetermined);
            newUndetermined.remove(next);

            // determine the actual costs
            List<Double> actualCosts = new ArrayList<>(newRotamers.size());
            for (Rotamer r : newRotamers)
                actualCosts.add( actualCost(r, this) );

            // calculate the components of the heuristic function that will be common to all the children
            // commonMap maps an undetermined rotamer r at position i (call it i_r) to:
            // E_self(i_r) +
            //
            // sum over j, where j is a determined position, of:
            //   E_interaction(i_r, j_s)
            //    --> normally, j would run over all determined positions; here, we only go over all common determined positions
            //
            // sum over k, where k != i, of:
            //   min over u, where u is a rotamer at k (call it k_u)
            //     E_interaction(i_r,k_u) --> excludes incompatibles
            Map<Rotamer,Double> commonMap = new HashMap<>();
            for (Integer i : newUndetermined)
                {
                    r:
                    for (Rotamer r : rotamerSpace.get(i))
                        {
                            // self-energy
                            double energy = rotamerSelfEnergies.get(r);
                            
                            // interaction energies with common determined rotamers
                            for (Rotamer s : rotamers) {
                                Double e = rotamerInteractionEnergies.get(new RotamerPair(r,s));
                                if ( e == null ) continue r; // don't include incompatible pairs
                                energy += e; 
                            }

                            // interaction energies with undetermined rotamers
                            double[] minima = energyMinima.get(r);
                            for (Integer k : newUndetermined)
                                if ( k != i )
                                    energy += minima[k];

                            // add this component to the map
                            commonMap.put(r, energy);
                        }
                }

            // calculate the heuristic costs
            List<Double> heuristicCosts = new ArrayList<>(newRotamers.size());
            for (Rotamer newRotamer : newRotamers)
                {
                    double thisHeuristicCost = 0.0;
                    for (Integer i : newUndetermined)
                        {
                            double min = 0.0;
                            boolean first = true;
                            r:
                            for (Rotamer r : rotamerSpace.get(i))
                                {
                                    Double energy = commonMap.get(r);
                                    if ( energy == null )
                                        continue r;

                                    // add the extra interaction component for this child node
                                    Double e = rotamerInteractionEnergies.get(new RotamerPair(r,newRotamer));
                                    if ( e == null )
                                        continue r;
                                    energy += e;

                                    // replace the minimum if necessary
                                    if ( first )
                                        {
                                            min = energy;
                                            first = false;
                                        }
                                    else if ( energy < min )
                                        min = energy;
                                }
                            thisHeuristicCost += min;
                        }
                    heuristicCosts.add(thisHeuristicCost);
                }

            // create the new nodes
            for (int i=0; i < newRotamers.size(); i++)
                {
                    // get fields
                    List<Rotamer> thisRotamers = newRotamerLists.get(i);
                    // newUndetermined is the undetermined for all children and they will all share the same object
                    double thisActualCost = actualCosts.get(i);
                    double thisTotalCost = thisActualCost + heuristicCosts.get(i);
                    
                    // create the new field and add it to the output
                    Node newNode = new Node(thisRotamers, newUndetermined, thisActualCost, thisTotalCost);
                    // double check the heuristic cost
                    //System.out.println( newNode.heuristicCost() + " : " + heuristicCosts.get(i) + " : " + (newNode.heuristicCost()-heuristicCosts.get(i)) );
                    output.add(newNode);
                }

            return output;
        }

        /**
         * Compares two Nodes based on their function values.  
         * @param other The other Node.  
         */
        public int compareTo(Node other)
        {
            return Double.compare(this.totalCost, other.totalCost); 
        }

        @Override
        public String toString()
        {
            String rotamerString = "";
            for (int i=0; i < rotamerSpace.size(); i++)
                {
                    boolean found = false;
                    for (Rotamer r : rotamers)
                        {
                            if ( r.sequenceIndex == i )
                                {
                                    rotamerString += String.format("%5s ", r.description.split("_")[0].substring(0,5));
                                    found = true;
                                    break;
                                }
                        }
                    if ( !found )
                        rotamerString += String.format("%5s ", "---");
                }
            return String.format("%10.2f : %s", totalCost, rotamerString);
        }

        @Override
        public boolean equals(Object obj)
        {
            if ( obj == null )
                return false;
            if ( obj == this )
                return true;
            if ( !(obj instanceof Node) )
                return false;

            Node n = (Node)obj;
            // note: this is an order-sensitive comparison, but should be safe since we always expand in the same order
            if ( Objects.equals(rotamers, n.rotamers) )
                return true;
            return false;
        }

        @Override
        public int hashCode()
        {
            return Objects.hash(rotamers);
        }
    } // end of class Node
    
    /** This work unit will do parallel A* expansions around nodes. */
    public class SearchUnit implements WorkUnit
    {
        public final PriorityQueue<Node> queue;
        public final CountDownLatch latch;
        public final List<Node> batch;

        public SearchUnit(PriorityQueue<Node> queue, CountDownLatch latch, List<Node> batch)
        {
            this.queue = queue;
            this.latch = latch;
            this.batch = batch;
        }

        public Result call()
        {
            // enqueue the nodes
            for (Node n : batch)
                queue.add(n);

            if ( queue.size() > 0 )
                {
                    // pop the highest priority node
                    //System.out.println("expanding");
                    Node highestPriorityNode = queue.remove();
            
                    // expand the node
                    List<Node> children = highestPriorityNode.expand();

                    // mark solutions if any
                    List<Node> theseSolutions = new ArrayList<>();
                    for (Node n : children)
                        {
                            if ( n.isTarget() )
                                theseSolutions.add(n);
                        }
                    children.removeAll(theseSolutions);
                    solutions.addAll(theseSolutions);
                    //System.out.printf("%d children and %d solutions\n", children.size(), theseSolutions.size());

                    // put explored nodes in a set common to all threads (newNodes)
                    newNodes.addAll(children);
                }

            // count down the latch
            latch.countDown();

            // dummy result
            return null;
        }
    }
}
