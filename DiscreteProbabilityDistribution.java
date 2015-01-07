
import java.util.*;
import java.util.concurrent.*;
import com.google.common.collect.*;

/**
 * A utility class for drawing random numbers from an arbitrary
 * discrete probability distribution using the alias method.  Please
 * see <a href="https://en.wikipedia.org/wiki/Alias_method">Wikipedia</a>
 * for more details.  This is a generic class that is parameterized by
 * the type of outcome that the distribution is meant to represent.  For
 * example, the outcome might be a list of torsion angles.  This class is
 * immutable.
 *
 * Algorithm copied from <a href="http://www.keithschwarz.com/darts-dice-coins/">here</a>.
 *
 */
public class DiscreteProbabilityDistribution<E>
{
    /**
     * Contains the discrete probability distribution itself.  Maps each
     * outcome to a probability in the range 0.0 - 1.0.
     */
    private final double[] probability;
    
    /** The outcomes, mapped 1:1 with the probability array. */
    public final ImmutableList<E> outcomes;

    /** The alias table. */
    private final int[] alias;

    /**
     * The initial probabilities.
     */
    public final ImmutableList<Double> inputProbabilities;

    /**
     * Populates the distribution and verifies invariants.  Lists must be
     * the same size and the sum of the probabilities must equal 1.
     * @param outcomes possible outcomes
     * @param probabilities probability of each outcome
    */
    public DiscreteProbabilityDistribution(List<E> outcomes, List<Double> probabilities)
    {
        // check for nulls
        if (outcomes == null || probabilities == null || outcomes.size() == 0 || probabilities.size() == 0)
            throw new NullPointerException("Cannot create a DiscreteProbabilityDistribution from a null or empty List!");

        // check lists are equal sizes
        if ( outcomes.size() != probabilities.size() )
            throw new IllegalArgumentException("Must have the same number of outcomes and probabilities!");

        // check for all positive probabilities
        double sum = 0.0;
        for (Double d : probabilities)
            {
                if ( d < 0 )
                    throw new IllegalArgumentException("Cannot have a negative probability!");
                sum += d;
            }
        
        // initialize arrays
        probability = new double[probabilities.size()];
        alias = new int[probabilities.size()];
        this.outcomes = ImmutableList.copyOf(outcomes);
        inputProbabilities = ImmutableList.copyOf(probabilities);

        // normalize probabilities and make a copy of the probabilities list, since we will be changing it
        probabilities = new ArrayList<Double>(probabilities);
        for (int i=0; i < probabilities.size(); i++)
            probabilities.set(i, probabilities.get(i) / sum);

        // calculate average probability
        final double average = 1.0 / probabilities.size();

        // create work lists
        Deque<Integer> small = new ArrayDeque<>();
        Deque<Integer> large = new ArrayDeque<>();

        // populate the input stacks with the input probabilities
        for (int i=0; i < probabilities.size(); ++i)
            {
                // if the probability is below avearge add it to small list
                // otherwise add it to the large list
                if ( probabilities.get(i) >= average )
                    large.add(i);
                else
                    small.add(i);
            }

        /* As a note: in the mathematical specification of the algorithm, we
         * will always exhaust the small list before the big list.  However,
         * due to floating point inaccuracies, this is not necessarily true.
         * Consequently, this inner loop (which tries to pair small and large
         * elements) will have to check that both lists aren't empty.
         */
        while ( !small.isEmpty() && !large.isEmpty() )
            {
                /* Get the index of the small and the large probabilities. */
                int less = small.removeLast();
                int more = large.removeLast();

                /* These probabilities have not yet been scaled up to be such that
                 * 1/n is given weight 1.0.  We do this here instead.
                 */
                probability[less] = probabilities.get(less) * probabilities.size();
                alias[less] = more;

                /* Decrease the probability of the larger one by the appropriate
                 * amount.
                 */
                probabilities.set( more, (probabilities.get(more) + probabilities.get(less)) - average );

                /* If the new probability is less than the average, add it into the
                 * small list; otherwise add it to the large list.
                 */
                if (probabilities.get(more) >= 1.0 / probabilities.size())
                    large.add(more);
                else
                    small.add(more);
            }

        /* At this point, everything is in one list, which means that the
         * remaining probabilities should all be 1/n.  Based on this, set them
         * appropriately.  Due to numerical issues, we can't be sure which
         * stack will hold the entries, so we empty both.
         */
        while (!small.isEmpty())
            probability[small.removeLast()] = 1.0;
        while (!large.isEmpty())
            probability[large.removeLast()] = 1.0;
    }
    /**
     * Draws a random outcome based on the given distribution using the alias method.
     * Uses thread-local random numbers.
     * @return E the outcome of the random weighted draw 
     */
    public E getRandom()
    {
        // draw a thread-safe random number
        int column = ThreadLocalRandom.current().nextInt(probability.length);

        // generate biased coin toss
        boolean coinToss = ThreadLocalRandom.current().nextDouble() < probability[column];

        // based on the outcome, get the colum or its alias
        int result = coinToss ? column : alias[column];

        // return the corresponding result object
        return outcomes.get(result);
    }
    
    /**
     * Returns a short description of this distribution:
     * [outcome, adjusted probability]
     * Note that this is not the actual or normalized probability!
     * To see that, you will have to print out the probabilities array in the
     * constructor.  (We throw that way afterwards.)
     * EDIT: I added back the original probabilities.  We should throw this
     * away in production to save memory.
     * @return the list [outcome, adjusted probability]
     */
    @Override
    public String toString()
    {
        String returnString = "";
        for (int i=0; i < probability.length; i++)
            {
                E outcome = outcomes.get(i);
                double prob = inputProbabilities.get(i);
                returnString = returnString + String.format("{outcome: %s, prob%% = %.2f}", outcome.toString(), 100.0 * prob);
                if ( i < probability.length - 1 )
                    returnString = returnString + "\n";
            }
        return returnString;
    }

    /**
     * Returns up to maxNumberOfLines probabilities above some threshold for debugging.
     * The biggest probabilities are returned first.
     *
     * @param threshold the minimum probability to print
     * @param maxNumberOfLines how many lines to print
     * @return the brief debug string
     */
    public String toDebugString(double threshold, int maxNumberOfLines)
    {
        TreeMap<Double,String> map = new TreeMap<>(Collections.reverseOrder());
        for (int i=0; i < probability.length; i++)
            {
                E outcome = outcomes.get(i);
                double prob = inputProbabilities.get(i);
                if ( prob < threshold )
                    continue;
                map.put(100.0 * prob, outcome.toString());
            }

        int count = 0;
        String returnString = "[ ";
        for (Double prob : map.keySet())
            {
                count++;
                if ( count > maxNumberOfLines )
                    break;
                String outcome = map.get(prob);
                returnString = returnString + String.format("{outcome: %s, prob%% = %.6f}\n", outcome, prob);
            }
        return returnString + " ]";
    }

    /**
     * Returns the hash code for this distribution.
     * @return the hash code for this distribution
     */
    @Override
    public int hashCode()
    {
        return Objects.hash(probability,outcomes,alias,inputProbabilities);
    }

    // equals
    /**
     * Returns true if this DiscreteProbabilityDistribution is identical to another
     * distribution.  Uses Arrays.equals().
     * @return true if equal
     */
    @Override
    public boolean equals(Object obj)
    {
        if ( obj == null )
            return false;
        if ( obj == this )
            return true;
        if ( !(obj instanceof DiscreteProbabilityDistribution) )
            return false;

        DiscreteProbabilityDistribution<?> dist = (DiscreteProbabilityDistribution<?>)obj;
        if ( ! Arrays.equals(dist.probability, probability) ||
             ! dist.outcomes.equals(outcomes) ||
             ! Arrays.equals(dist.alias, alias) ||
             ! dist.inputProbabilities.equals(inputProbabilities) )
            return false;
        return true;
    }

    public int getSize() 
    {
	return outcomes.size();
    }

    /**
     * Tester class.
     */
    public static void main(String[] args)
    {
        // generate test distribution
        List<String> outcomes = new LinkedList<String>();
        outcomes.add("A");
        outcomes.add("B");
        outcomes.add("C");
        outcomes.add("D");
        outcomes.add("E");

        List<Double> probabilities = new LinkedList<Double>();
        probabilities.add(2.0);
        probabilities.add(1.0);
        probabilities.add(1.0);
        probabilities.add(1.0);
        probabilities.add(1.0);

        DiscreteProbabilityDistribution<String> dist = new DiscreteProbabilityDistribution<>(outcomes,probabilities);
        LinkedHashMap<String,Integer> results = new LinkedHashMap<>();
        for (String s : outcomes)
            results.put(s, 0);
        System.out.println(dist);
        System.out.println("Rolling...");
        for (int i=0; i < 600000; i++)
            {
                String thisOutcome = dist.getRandom();
                Integer numberOfHits = results.get(thisOutcome);
                numberOfHits = numberOfHits + 1;
                results.put(thisOutcome, numberOfHits);
            }
        System.out.println(results);
    }
}
