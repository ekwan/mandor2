import java.io.*;
import java.util.*;
import com.google.common.collect.*;

/**
 * A class that is used to parse the output from a Tinker call to analyze.
 * It creates an energy breakdown of the total energy of a peptide into residue
 * components by treating the energy of a component as a self residue and an
 * interaction term.  Residue-residue interactions are divided equally between
 * the two residues.  There is no "backbone"--the backbone portion of each residue
 * is considered to be part of the residue.
 */
public class TinkerAnalyzeOutputFile extends OutputFileFormat
{
    /** the energy of each residue in kcal indexed by residue number */
    public final List<Double> energyByResidue;

    /** the total energy of the peptide in kcal */
    public final double totalEnergy;
    
    /**
     * Analyzes the residue-by-residue energy in a peptide and stores the result.
     * @param filename the name of the file containing the output from a tinker call to analyze (option D, "details")
     * @param peptide the corresponding peptide
     */
    public TinkerAnalyzeOutputFile(String filename, Peptide peptide)
    {
        // read file
        super(filename);

        // initialize array
        double[][] energyByResidue = new double[peptide.sequence.size()][peptide.sequence.size()];

        // make a map from atom numbers to residue
        List<Residue> sequence = peptide.sequence;
        Map<Integer,Residue> residueMap = new HashMap<>();
        
        for (Residue r : sequence)
        {
            for (Atom a : r.atoms) 
            {
                int atomNumber = peptide.getAtomNumber(a);
                if ( atomNumber <= 0 )
                    throw new IllegalArgumentException("atom not found");
                else if ( residueMap.containsKey(atomNumber) )
                    throw new IllegalArgumentException("duplicate atom during tinker analysis");
                residueMap.put(atomNumber, r);
            }
        }

        // sum of all individual energies used to check if method works
        Double totalEnergy = 0.0;

        // parse all interactions
        for (List<String> line : fileContents) 
            {
                String heading = null;
                if ( line.size() > 0 )
                    heading = line.get(0);
		        String heading2 = null;
                if ( line.size() > 1 ) 
                    heading2 = line.get(1);
                if ( ( heading.equals("Bond") && heading2.equals("Stretching")  ) ||
                     ( heading.equals("Angle") && heading2.equals("Bending")    ) ||
                     ( heading.equals("Improper") && heading2.equals("Torsion") )    )
                     continue;
                else if ( heading.equals("Bond") || heading.equals("PiTors") || heading.equals("VDW-Hal") || heading.equals("VDW-LJ") || heading.equals("Charge") )
                    {
			            int atomNumber1 = getInt(line.get(1));
                        int atomNumber2 = getInt(line.get(2));
                        
                        Pair<Integer, Integer> interactionClassification = classifyInteraction(sequence, residueMap, atomNumber1, atomNumber2);
                        
                        int residue1 = interactionClassification.getFirst();
                        int residue2 = interactionClassification.getSecond();
                        if ( heading.equals("PiTors") )
                            energyByResidue[residue1][residue2] += Double.parseDouble(line.get(4));
                        else if ( heading.equals("Charge") )
                            energyByResidue[residue1][residue2] += Double.parseDouble(line.get(6));
                        else
                            energyByResidue[residue1][residue2] += Double.parseDouble(line.get(5));
                    }
                else if (heading.equals("Angle") || heading.equals("Angle-IP") || heading.equals("StrBend"))
                    {
                        int atomNumber1 = getInt(line.get(1));
                        int atomNumber2 = getInt(line.get(2));
                        int atomNumber3 = getInt(line.get(3));

                        Pair<Integer,Integer> interactionClassification = classifyInteraction(sequence, residueMap, atomNumber1, atomNumber2, atomNumber3);
                        int residue1 = interactionClassification.getFirst();
                        int residue2 = interactionClassification.getSecond();

                        energyByResidue[residue1][residue2] = energyByResidue[residue1][residue2] + Double.parseDouble(line.get(6));
                    }
                else if (heading.equals("Torsion") || heading.equals("O-P-Bend") || heading.equals("Improper") )
		            {
                        int atomNumber1 = getInt(line.get(1));
                        int atomNumber2 = getInt(line.get(2));
                        int atomNumber3 = getInt(line.get(3));
                        int atomNumber4 = getInt(line.get(4));

                        Pair<Integer,Integer> interactionClassification = classifyInteraction(sequence, residueMap, atomNumber1, atomNumber2, atomNumber3, atomNumber4);
                        int residue1 = interactionClassification.getFirst();
                        int residue2 = interactionClassification.getSecond();
                        energyByResidue[residue1][residue2] = energyByResidue[residue1][residue2] + Double.parseDouble(line.get(line.size()-1));
		            }
                else if (heading.equals("TorTor"))
                    {
                        int atomNumber1 = getInt(line.get(1));
                        int atomNumber2 = getInt(line.get(2));
                        int atomNumber3 = getInt(line.get(3));
                        int atomNumber4 = getInt(line.get(4));
                        int atomNumber5 = getInt(line.get(5));

                        Pair<Integer,Integer> interactionClassification = classifyInteraction(sequence, residueMap, atomNumber1, atomNumber2, atomNumber3, atomNumber4, atomNumber5);
                        int residue1 = interactionClassification.getFirst();
                        int residue2 = interactionClassification.getSecond();
                        energyByResidue[residue1][residue2] = energyByResidue[residue1][residue2] + Double.parseDouble(line.get(line.size()-1));
                    }
	    	    else if (heading.equals("M-Pole"))
                    {
                        int atomNumber1 = getInt(line.get(1));
                        int atomNumber2 = getInt(line.get(2));
                        Pair<Integer, Integer> interactionClassification = classifyInteraction(sequence,residueMap, atomNumber1, atomNumber2);
                        int residue1 = interactionClassification.getFirst();
                        int residue2 = interactionClassification.getSecond();
                        energyByResidue[residue1][residue2] = energyByResidue[residue1][residue2] + Double.parseDouble(line.get(4));
			            energyByResidue[residue1][residue2] = energyByResidue[residue1][residue2] + Double.parseDouble(line.get(5));
                    }
		        else if ( heading.equals("Solvate"))
                    {
                        int atomNumber = getInt(line.get(1));
                        Pair<Integer, Integer> interactionClassification = classifyInteraction(sequence,residueMap, atomNumber);
                        int residue1 = interactionClassification.getFirst();
                        int residue2 = interactionClassification.getSecond();
                        if ( residue1 != residue2 )
                            throw new IllegalArgumentException("mismatch");
                        energyByResidue[residue1][residue2] = energyByResidue[residue1][residue2] + Double.parseDouble(line.get(4));
                    }
                else if (heading.equals("Total") && line.get(1).equals("Potential") && line.get(2).equals("Energy"))
                    {
                        totalEnergy = Double.parseDouble(line.get(4));
                        //System.out.println("tinker energy: " + totalEnergy);
                    }
            }

        // energyOfResidue contains the energy of each residue
        List<Double> energyOfResidue = new ArrayList<>();
	    for (int i=0; i<sequence.size(); i++)
	        energyOfResidue.add(0.0);
	
	    // check point energy
	    //Double checkEnergy = 0.0;

        for (int i = 0; i < sequence.size() ; i++)
            for (int j = 0; j < sequence.size(); j++)
            {
		        //checkEnergy = checkEnergy + energyByResidue[i][j];

                if (i == j)
                {
                    // self energy term
		            Double currentEnergyTerm = energyOfResidue.get(i);
                    energyOfResidue.set(i,currentEnergyTerm + energyByResidue[i][j]);
                }
                else
                {
                    // represents interaction terms
		            Double currentEnergyTermI = energyOfResidue.get(i);
		            Double currentEnergyTermJ = energyOfResidue.get(j);
                    energyOfResidue.set(i,currentEnergyTermI + energyByResidue[i][j] / 2.0);
                    energyOfResidue.set(j,currentEnergyTermJ + energyByResidue[i][j] / 2.0);
                }
            }

	    // energy total from individual resiudes used to check
        Double totalEnergyFromResidues = 0.0;
        for (Double e : energyOfResidue)
	        totalEnergyFromResidues = totalEnergyFromResidues + e;
        
        if (Math.abs(totalEnergyFromResidues-totalEnergy)>0.15)
            throw new IllegalArgumentException("total energy does not equal individual energies. \nEnergy of residues is " + totalEnergyFromResidues + "\nEnergy from file is " + totalEnergy + "\ndifference: " + (totalEnergyFromResidues - totalEnergy));

        //System.out.println("check energy:               " + checkEnergy);
	    //System.out.println("total energy from residues: " + totalEnergyFromResidues);
	    //System.out.println("tinker energy:              " + totalEnergy);

        this.totalEnergy = totalEnergy;
        this.energyByResidue = ImmutableList.copyOf(energyOfResidue);
    }
	
    // private method that classifies a set of atom numbers into the corresponding residues (returned as index in sequence list)
    private Pair<Integer,Integer> classifyInteraction(List<Residue> sequence, Map<Integer,Residue> residueMap, int ... atomNumbers)
    {
        Residue residue1 = residueMap.get(atomNumbers[0]);
        if ( residue1 == null )
            throw new NullPointerException("can't have a null residue");
            
        Residue residue2 = null;

        for (int i=1; i < atomNumbers.length; i++)
            {
                if ( ! residueMap.get(atomNumbers[i]).equals(residue1) )
                    {
                        residue2 = residueMap.get(atomNumbers[i]);
                        break;
                    }
            }
        
        if ( residue2==null )
            return new Pair<Integer,Integer>(sequence.indexOf(residue1), sequence.indexOf(residue1));
        else
            return new Pair<Integer,Integer>(sequence.indexOf(residue1), sequence.indexOf(residue2));
    }
    
    /**
     * Returns the int corresponding to a given string.  Strips all non-numeric characters from the string.
     * @param s the the string to parse
     * @return the int corresponding to s
     */
    private int getInt(String s)
    {
        String newString = null;
        try { newString = s.replaceAll("[^\\d]",""); }
        catch (Exception e) { System.out.println(s); throw e; }
        return Integer.parseInt(newString);
    }
 
    @Override
    public String toString()
    {
        return energyByResidue.toString() + "\nTotal Energy:  " + totalEnergy;
    }

    @Override
    public int hashCode()
    {
        return Objects.hash(stringRepresentation, fileContents, energyByResidue, totalEnergy);
    }

    @Override
    public boolean equals(Object obj)
    {
        if ( obj == null )
            return false;
        if ( obj == this )
            return true;
        if ( !(obj instanceof TinkerAnalyzeOutputFile) )
            return false;

        TinkerAnalyzeOutputFile o = (TinkerAnalyzeOutputFile)obj;
        if ( stringRepresentation.equals(o.stringRepresentation) &&
             fileContents.equals(o.fileContents) &&
             energyByResidue.equals(o.energyByResidue) &&
             totalEnergy == o.totalEnergy )
            return true;
        return false;
    }
}
