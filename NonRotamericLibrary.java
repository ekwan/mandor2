import java.util.*;
import java.io.*;
import com.google.common.collect.*;

/**
 * A class used to read in rotamer data from the <a href="http://dunbrack.fccc.edu/bbdep2010/">Dunbrack library</a>
 * It allows you to access the DiscreteProbabilityDistribution associated with a set of backbone angles, phi and psi.
 * This class is immutable.  Rotamers below the threshold value Settings.ROTAMER_LIBRARY_THRESHOLD will
 * be ignored.  We assume that phi and psi span the interval [-180.0, 180.0] with step 10.0. <p>
 * This library is for non-rotameric sidechains (i.e., those with terminal sp2-sp3 torsions).
 */
public class NonRotamericLibrary extends RotamerLibrary
{
    /**
     * Utility class for the probabilistic outcomes.  Contains a list of standard torsion angles
     * and a distribution for the last non-rotameric angle.  This class is immutable.
     */
    public class NonRotamericAngles
    {
        /** X1, X2, ..., Xn-1 torsion angles in degrees */
        private final List<Double> standardTorsionAngles;

        /** a discrete distribution for Xn */
        private final DiscreteProbabilityDistribution<Double> nonRotamericTorsionAngleDistribution;

        /** constructs an immutable instance of this class */
        public NonRotamericAngles(List<Double> chis, DiscreteProbabilityDistribution<Double> dpd)
        {
            standardTorsionAngles = ImmutableList.copyOf(chis);
            nonRotamericTorsionAngleDistribution = dpd;
        }

        /**
         * returns X1, X2, ..., Xn-1 torsion angles in degrees 
         * @return X1, X2, ..., Xn-1 torsion angles in degrees
         */
	    public List<Double> getRotamericAngles()
        {
	        return standardTorsionAngles;
	    }
	
        /**
         * returns the discrete probability distribution for the non-rotameric torsion angle
         * @return the discrete probability distribution for the non-rotameric torsion angle (degrees)
         */
	    public DiscreteProbabilityDistribution<Double> getDPD()
        {
            return nonRotamericTorsionAngleDistribution;
	    }

        @Override
        public String toString()
        {
            return String.format("X1...Xn-1: %s\nXn: %s\n", standardTorsionAngles.toString(), nonRotamericTorsionAngleDistribution.toString());
        }

        @Override
        public int hashCode()
        {
            return Objects.hash(standardTorsionAngles, nonRotamericTorsionAngleDistribution);
        }

        @Override
        public boolean equals(Object obj)
        {
            if ( obj == null )
                return false;
            if ( ! (obj instanceof NonRotamericAngles) )
                return false;
            final NonRotamericAngles other = (NonRotamericAngles)obj;
            return (Objects.equals(standardTorsionAngles, other.standardTorsionAngles) &&
                    Objects.equals(nonRotamericTorsionAngleDistribution, other.nonRotamericTorsionAngleDistribution));
        }
    }
				
    /**
      * Data is stored in this nested structure:
      * BackboneAngles(phi,psi) --> DiscreteProbabilityDistribution( NonRotamericAngles, Double probabilities )
      * where NonRotamericAngles is a List<Double> of chis and a DiscreteProbabilityDistribution
      * in the last non-rotameric angle.
      */
    private final Map<RotamerLibrary.Angles, DiscreteProbabilityDistribution<NonRotamericAngles>> dataset; 
    
    /** The number n-1 of normal chi angles: X1, X2, ..., Xn-1 */
    private final Integer numberOfSidechainTorsions;

    /** The type of amino acid these data are for. */
    private final AminoAcid aminoAcid;

    /**
     * creates a RotamericLibrary by reading in filename associated with amino acid
     * @param aminoAcid the amino acid for which to create rotameric library
     */
    public NonRotamericLibrary(AminoAcid aminoAcid, String filename)
    {
	    // read in library for specific amino acid
	    // creates LinkedHashMap between backbone angles and discrete probability distribution
        this.aminoAcid = aminoAcid;

        // read in entire file
        Map<RotamerLibrary.Angles, DiscreteProbabilityDistribution<NonRotamericAngles>> tempDataset = new LinkedHashMap<>();
        Integer tempNumberOfSidechainTorsions = null;
        try (BufferedReader thisFile = new BufferedReader(new FileReader(filename)))
            {
                // make sure this is the right kind of file
                if ( filename.indexOf("rotamers") > -1 )
                    throw new IllegalArgumentException("Should not be using RotamericLibrary to read standard rotameric side chains!");
		        
                // stores which fields have regular chi data in them
                List<Integer> chiFieldNumbers = new LinkedList<>();
                int probabilityFieldIndex = 0;

                // stores which fields have non-rotameric chi data in them
                List<Double> nonRotamericChiLabels = new LinkedList<>();
                List<Integer> nonRotamericChiFieldNumbers = new LinkedList<>();

                Double lastPhi = null;
                Double lastPsi = null;

                List<NonRotamericAngles> tempNRA = new LinkedList<>();   // the domain for this line
                List<Double> tempProbabilities   = new LinkedList<>();   // the range (i.e., probability) for this line of (X1, X2, ..., Xn-1, non-rotameric distribution)

                String currentLine = null;
                while ( (currentLine = thisFile.readLine()) != null )
                    {
			            // parse currentLine to access data fields
			            String[] parts = currentLine.split("\\s+");
		    
			            // valid data is stored on lines that contain the
                        // abbreviated amino acid name as the first field
                        if ( currentLine.indexOf("Number of chi angles treated as discrete") > -1 )
                            {
                                // records the number of torsions we need to store (n-1)
                                tempNumberOfSidechainTorsions = Integer.valueOf(parts[parts.length-1]);
                                continue;
                            }
                        else if ( currentLine.indexOf("chi1Val") > -1 )
                            {
                                // compute which field have chis(1,2,...,n-1)
                                int currentField = 4;
                                for (int i=currentField; i <= 4+tempNumberOfSidechainTorsions; i++)
                                    currentField++;
                                //System.out.println("chi val start:" + currentField);
                                for (int i=currentField; i < currentField+tempNumberOfSidechainTorsions; i++)
                                    chiFieldNumbers.add(i);
                                //System.out.println("Chi field numbers: " + chiFieldNumbers.toString());

                                // compute which fields have chi(n) labels
                                currentField = chiFieldNumbers.get(0) - 1;
                                probabilityFieldIndex = currentField;
                                //System.out.println("prob:" + currentField);
                                currentField = currentField + 2*chiFieldNumbers.size() + 1;
                                //System.out.println("start:" + currentField);
                                for (int i=currentField; i < parts.length-1; i++)
                                    {
                                        nonRotamericChiLabels.add(Double.valueOf(parts[i+1]));
                                        nonRotamericChiFieldNumbers.add(i);
                                    }

                                //System.out.println("non-rotameric fields: " + nonRotamericChiFieldNumbers.toString());
                                //System.out.println("labels: " + nonRotamericChiLabels);
                                //System.out.println( nonRotamericChiFieldNumbers.size() + " / " + nonRotamericChiLabels.size() );
                                continue;
                            }
			            else if (!parts[0].equals(aminoAcid.toString().toUpperCase()))
                            continue;

                        // ignore low-probability rotamers
                        Double probability = Double.valueOf(parts[probabilityFieldIndex]);
			            if ( probability < Settings.ROTAMER_LIBRARY_THRESHOLD )
                            continue;

			            // read backbone angles with parts.get(1) and parts.get(2)
			            Double currPhi = Double.parseDouble(parts[1]);
			            Double currPsi = Double.parseDouble(parts[2]);
                        if ( lastPhi == null )
                            lastPhi = currPhi;
                        if ( lastPsi == null )
                            lastPsi = currPsi;

                        // if we are working on a new phi,psi pair, then
                        // load the last block into the database
                        // this assumes entries are consecutive for one pair of (phi,psi)
                        if ( !currPhi.equals(lastPhi) || !currPsi.equals(lastPsi))
                            {
                                RotamerLibrary.Angles bba = new RotamerLibrary.Angles(lastPhi,lastPsi);
                                
                                // create the distribution that stores entries for this entire block
                                DiscreteProbabilityDistribution<NonRotamericAngles> outerDPD = new DiscreteProbabilityDistribution<>(tempNRA,tempProbabilities); 
                                tempDataset.put(bba, outerDPD);
                                tempNRA = new LinkedList<>();
                                tempProbabilities = new LinkedList<>();
                            }

                        // read chi angles
                        List<Double> regularChis = new LinkedList<>();
                        for (Integer i : chiFieldNumbers)
                            regularChis.add(Double.valueOf(parts[i]));

                        // read non-rotameric chi distribution
                        List<Double> temp = new LinkedList<>();
                        for (Integer i : nonRotamericChiFieldNumbers)
                            temp.add(Double.valueOf(parts[i]));

                        // create DiscreteProbabilityDistribution that stores that non-rotameric distribution for this line
                        // note: everything points to the same list of labels!  efficiency here we come!
                        DiscreteProbabilityDistribution<Double> innerDPD = new DiscreteProbabilityDistribution<>(nonRotamericChiLabels, temp); 

                        // create NonRotamericAngles
                        NonRotamericAngles thisNRA = new NonRotamericAngles(regularChis,innerDPD);

                        // add this line to the temporary list
                        tempNRA.add(thisNRA);
                        tempProbabilities.add(probability);

                        // reset for next round
                        lastPhi = currPhi;
                        lastPsi = currPsi;
		        }

                // edge case: store last block of data
		        RotamerLibrary.Angles bba = new RotamerLibrary.Angles(lastPhi,lastPsi);

		        // create the distribution that stores entries for this entire block
		        DiscreteProbabilityDistribution<NonRotamericAngles> outerDPD = new DiscreteProbabilityDistribution<>(tempNRA,tempProbabilities);
		        tempDataset.put(bba, outerDPD);
            }
        catch (IOException e)
            {
                e.printStackTrace();
                System.exit(1);
            }
            
        // make defensive and immutable copy
        dataset = ImmutableMap.copyOf(tempDataset);
        numberOfSidechainTorsions = tempNumberOfSidechainTorsions;
    }

    /**
     * returns the DiscreteProbabilityDistribution associated with phi and psi
     * rounds psi and phi to nearest multiple of 10
     * @param phi a Double with the corresponding phi from the backbone
     * @param psi a Double with the corresponding psi from the backbone
     * @return a DPD object with list of chis as outcomes and associated probabilities for a phi and psi
     */
    public DiscreteProbabilityDistribution<NonRotamericAngles> get(Double phi, Double psi)
    {
        if (phi > 180.0 || phi < -180.0 || psi > 180.0 || psi < -180.0)
            throw new IllegalArgumentException("psi and phi must be between -180 and 180");
        
        // round to nearest 10 and return the appropriate data
        double phi_rounded = roundTo10(phi);
        double psi_rounded = roundTo10(psi);
        
        DiscreteProbabilityDistribution<NonRotamericAngles> dpd = dataset.get(new RotamerLibrary.Angles(phi_rounded,psi_rounded));
        
        if ( dpd == null )
	    {
		    System.out.println("phi_rounded is " + phi_rounded);
		    System.out.println("psi_rounded is " + psi_rounded);
		    System.out.println(aminoAcid);
		    throw new NullPointerException("data not found!");
	    }
        return dpd;
    }

    /**
     * Rounds a double to the nearest 10.
     * @param number the number to be rouned
     * @return the double that is equal to the value of the nearest multiple of 10
     */
    public static double roundTo10(double number)
    {
        // prevents "negative zero"
        double rounded = Math.rint(number/10.0)*10.0 + 0.0;
        //System.out.println(number + " -> " + rounded);
        return rounded;
    }

    @Override
    public String toString()
    {
	    return String.format("NonRotamericLibrary for %s comprising %d entries.", aminoAcid.toString(), dataset.size());
    }
    
    @Override
    public boolean equals(Object obj)
    {
        if ( obj == null )
                return false;
        if ( obj == this )
            return true;
        if ( !(obj instanceof RotamericLibrary) )
            return false;

        NonRotamericLibrary l = (NonRotamericLibrary) obj;
        if ( dataset.equals(l.dataset) &&
             numberOfSidechainTorsions.equals(l.numberOfSidechainTorsions) &&
             aminoAcid == l.aminoAcid )
            return true;
        return false;
    }

    @Override
    public int hashCode()
    {
        return Objects.hash(dataset, aminoAcid, numberOfSidechainTorsions);
    }

    /** Tests the functionality of this library. */
    public static void main(String[] args)
    {
	    NonRotamericLibrary rotLib1 = new NonRotamericLibrary(AminoAcid.ASN, "rotamer_library/asn.bbdep.densities.lib");
	    System.out.println(rotLib1);
        System.out.println(rotLib1.get(81.0,60.1).toString());
    }
}
