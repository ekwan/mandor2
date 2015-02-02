import java.util.*;
import java.io.*;
import org.apache.commons.math3.geometry.euclidean.threed.*;
import org.jgrapht.*;
import org.jgrapht.graph.*;
import org.jgrapht.alg.*;
import com.google.common.collect.*;
import java.util.concurrent.*;

/**
 * This class collects together some methods for summarizing rotamer chi distributions.
 */
public class RotamerSummarizer
{
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

    /** This class is not instantiable. */
    private RotamerSummarizer()
    {
        throw new IllegalArgumentException("not instantiable");
    }

    /**
     * Given a residue, return the possible rotamers.  Rotamers are based on the current
     * phi and psi values.  Note that this might produce some rotamers that clash with
     * other parts of the peptide.
     *
     * Will throw an exception if there are no library data for the specified residue.
     *
     * Rotamers beneath PROBABILITY_THRESHOLD will be ignored.  Non-rotameric degrees of freedom
     * will be summarized into several angles using BIN_SIZE and BIN_THRESHOLD.
     *
     * @param residue the input residue
     * @return a nested list of all the angles chi1, chi2, ..., chiN
     */
    public static List<List<Double>> getPossibleRotamers(Residue residue)
    {
        AminoAcid aminoAcid = residue.aminoAcid;
        AminoAcid.RotamerType rotamerType = aminoAcid.rotamerType;
        double omega = residue.omega.getDihedralAngle();
        double phi = residue.phi.getDihedralAngle();
        double psi = residue.psi.getDihedralAngle();

        if ( rotamerType == AminoAcid.RotamerType.HAS_NO_ROTAMERS )
            return ImmutableList.of();
        else if ( rotamerType == AminoAcid.RotamerType.SPECIAL )
            throw new IllegalArgumentException("no data for special amino acids");
        if (rotamerType == AminoAcid.RotamerType.IS_ROTAMERIC)
            {
                RotamericLibrary rotLib = (RotamericLibrary)RotamerDatabase.getLibrary(aminoAcid, omega);
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
                if ( prunedRotamers.size() == 0 )
                    throw new IllegalArgumentException("pruned rotamers cannot be empty");
                return prunedRotamers;
            }

        else if (rotamerType == AminoAcid.RotamerType.NON_ROTAMERIC)
            {
                NonRotamericLibrary nRotLib = (NonRotamericLibrary)RotamerDatabase.getLibrary(aminoAcid, omega);
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
            }

        // should be unreachable
        throw new IllegalArgumentException("unreachable");
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

    /** For testing. */
    public static void main(String[] args)
    {
        DatabaseLoader.go();
        List<ProtoAminoAcid> sequence = ProtoAminoAcidDatabase.getSpecificSequence("arg","met","standard_ala","gly","d_proline", "gly", "phe", "val", "hd", "l_pro");
        Peptide peptide = PeptideFactory.createPeptide(sequence);
        
        for (int i=0; i < peptide.sequence.size(); i++)
            {
                peptide = BackboneMutator.mutateOmega(peptide, i);
                peptide = BackboneMutator.mutatePhiPsi(peptide, i);
                peptide = RotamerMutator.mutateChis(peptide, i);
            }
        peptide = PeptideFactory.setHairpinAngles(peptide);

        Residue residue = peptide.sequence.get(0);
        System.out.println(getPossibleRotamers(residue));
    }
}
