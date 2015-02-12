import java.io.*;
import java.util.*;
import com.google.common.collect.*;
import org.apache.commons.math3.geometry.euclidean.threed.*;
import org.jgrapht.*;
import org.jgrapht.graph.*;
import org.jgrapht.alg.*;

/**
 * This class provides static methods for calculating the interactions between
 * atoms in a molecule using the OPLS forcefield.
 */
public class OPLScalculator implements Singleton
{
    /**
     * the coulombic constant for calculating electrostatic interactions, which is calculated as
     * e*e/(4*pi*epsilon_0) in A and kcal/mol
     */
    //public static final double COULOMB_CONSTANT = 332.3814783;
    public static final double COULOMB_CONSTANT = 332.06;

    /** distances that are less than this value will be set to this value to avoid blowups */
    public static final double MIN_DISTANCE = 1.0;

    /** cannot be instantiated */
    private OPLScalculator()
    {
        throw new IllegalArgumentException("not instantiable");
    }

    /**
     * Gets the OPLS partial charge.
     * @param type the OPLS atom type
     * @return the partial charge in fractions of an electron charge
     */
    public static Double getCharge(int type)
    {
        Double charge = OPLSforcefield.CHARGE_MAP.get(type);
        if ( charge == null )
            throw new NullPointerException("charge not found for type " + type);
        return charge;
    }

    /**
     * Gets the VDW distance.
     * @param classNumber the OPLS atom class
     * @return the VDW potential distance in angstroms
     */
    public static Double getVDWDistance(Integer classNumber)
    {
        Double distance = OPLSforcefield.VDW_DISTANCE_MAP.get(classNumber);
        if ( distance == null )
            throw new NullPointerException("vdw not found for class " + classNumber);
        return distance;
    }

    /**
     * Gets the VDW distance.
     * @param classNumber the OPLS atom class
     * @return the VDW well depth in kcal/mol
     */
    public static Double getVDWDepth(Integer classNumber)
    {
        Double depth = OPLSforcefield.VDW_DEPTH_MAP.get(classNumber);
        if ( depth == null )
            throw new NullPointerException("vdw not found for class " + classNumber);
        return depth;
    }

    /**
     * Get the OPLS atom class.
     * @param type the OPLS atom type
     * @return the OPLS atom class
     */
    public static Integer getOPLSClass(int type)
    {
        Integer classNumber = OPLSforcefield.CLASS_MAP.get(type);
        if ( classNumber == null )
            throw new NullPointerException("class not found for type " + type);
        return classNumber;
    }

    /**
     * Computes the non-bonded OPLS interactions in a molecule between all atom pairs.
     * @param molecule the molecule to analyze
     * @return list of all non-bonded interactions
     */
    public static List<Interaction> getInteractions(Molecule molecule)
    {
        // get fields
        List<Atom> contents = molecule.contents;
        SimpleWeightedGraph<Atom,DefaultWeightedEdge> connectivity = molecule.connectivity;

        // generate list of valid atom pairs
        int estimatedSize = contents.size() * (contents.size()-1) / 2;
        List<Interaction> interactions = new ArrayList<Interaction>(estimatedSize);
        for (int i=0; i < contents.size(); i++)
            {
                Atom atom1 = contents.get(i);
                Integer atomType1 = atom1.type2;
                Integer atomClass1 = getOPLSClass(atomType1);

                double charge1 = getCharge(atomType1);
                double vdw_distance1 = getVDWDistance(atomClass1);
                double vdw_depth1 = getVDWDepth(atomClass1);

                for (int j=i+1; j < contents.size(); j++)
                    {
                        Atom atom2 = contents.get(j);
                        Integer atomType2 = atom2.type2;
                        Integer atomClass2 = getOPLSClass(atomType2);

                        double charge2 = getCharge(atomType2);
                        double vdw_distance2 = getVDWDistance(atomClass2);
                        double vdw_depth2 = getVDWDepth(atomClass2);

                        // calculate the graph-theoretic distance
                        DijkstraShortestPath<Atom,DefaultWeightedEdge> path = new DijkstraShortestPath<>(connectivity, atom1, atom2, 3.0);
                        List<DefaultWeightedEdge> pathEdges = path.getPathEdgeList();

                        // scale 1,4-interactions by 50%
                        double scaling = 1.0;
                        if ( pathEdges != null && pathEdges.size() < 3 )
                            continue;
                        else if ( pathEdges != null && pathEdges.size() == 3 )
                            scaling = 0.5;

                        // calculate the coulombic energy
                        double distance = Vector3D.distance(atom1.position, atom2.position);
                        //double vdwSum = 0.8 * (vdw_distance1 + vdw_distance2);
                        //if ( distance < vdwSum )
                        //    distance = vdwSum;
                        if ( distance < MIN_DISTANCE )
                            distance = MIN_DISTANCE;
                        double electrostatic = (charge1 * charge2 * COULOMB_CONSTANT) / (distance * distance);

                        // apply the combining rules for epsilon and sigma if necessary
                        double sigma = vdw_distance1;
                        double epsilon = vdw_depth1;

                        if ( atomClass1 != atomClass2 )
                            {
                                sigma = Math.sqrt(vdw_distance1 * vdw_distance2);
                                epsilon = Math.sqrt(vdw_depth1 * vdw_depth2);
                            }

                        // calculate the steric energy
                        double temp = Math.pow(sigma / distance, 6);
                        double steric = 4.0 * epsilon * temp * ( temp - 1.0 );

                        // apply scaling
                        double energy = ( electrostatic + steric ) * scaling;
                        
                        // for debugging
                        /*System.out.printf("Interaction for %s and %s (distance = %.4f):\n", molecule.getAtomString(atom1), molecule.getAtomString(atom2), distance);
                        System.out.printf("type = %3d  q = %7.4f  sigma = %7.4f  epsilon = %7.4f\n", atomType1, charge1, vdw_distance1, vdw_depth1);
                        System.out.printf("type = %3d  q = %7.4f  sigma = %7.4f  epsilon = %7.4f\n", atomType2, charge2, vdw_distance2, vdw_depth2);
                        System.out.printf("combined sigma = %7.4f  epsilon = %7.4f\n", sigma, epsilon);
                        System.out.printf("Electrostatic component: %.4f kcal\n", electrostatic);
                        System.out.printf("Steric component:        %.4f kcal\n", steric);
                        System.out.printf("Scaling factor:          %.4f\n", scaling);
                        System.out.printf("Total Potential Energy:  %.4f kcal\n", energy);
                        */

                        /*String description = String.format("distance = %5.2f    q1: %3d, %7.4f   q2 %3d, %7.4f   coulomb: %10.2f",
                                                            distance, atomType1.intValue(), charge1, atomType2.intValue(), charge2, electrostatic);

                        if ( scaling != 1.0 )
                            description += " [scaled]";
                        */

                        // create the interaction
                        Set<Atom> atomList = ImmutableSet.of(atom1, atom2);
                        //Interaction interaction = new Interaction(atomList, energy, description);
                        Interaction interaction = new Interaction(atomList, energy);
                        interactions.add(interaction);
                    }
            }
        
        // return the result
        return ImmutableList.copyOf(interactions);
    }

    /**
     * Computes the OPLS non-bonded interaction energy between two rotamers.
     * @param rotamer1 the first rotamer
     * @param rotamer2 the second rotamer
     * @return the interaction energy in kcal
     */
    public static double getInteractionEnergy(Rotamer rotamer1, Atom extraAtom1, Rotamer rotamer2, Atom extraAtom2)
    {
        double energy = 0.0;
        ArrayList<Atom> list1 = new ArrayList<Atom>(rotamer1.atoms);
        if ( extraAtom1 != null )
            list1.add(extraAtom1);
        ArrayList<Atom> list2 = new ArrayList<Atom>(rotamer2.atoms);
        if ( extraAtom2 != null )
            list2.add(extraAtom2);
        for (Atom atom1 : list1)
            {
                Integer atomType1 = atom1.type2;
                Integer atomClass1 = getOPLSClass(atomType1);

                double charge1 = getCharge(atomType1);
                double vdw_distance1 = getVDWDistance(atomClass1);
                double vdw_depth1 = getVDWDepth(atomClass1);
                
                for (Atom atom2 : list2)
                    {
                        Integer atomType2 = atom2.type2;
                        Integer atomClass2 = getOPLSClass(atomType2);

                        double charge2 = getCharge(atomType2);
                        double vdw_distance2 = getVDWDistance(atomClass2);
                        double vdw_depth2 = getVDWDepth(atomClass2);

                        // calculate the coulombic energy
                        double distance = Vector3D.distance(atom1.position, atom2.position);
                        //double vdwSum = 0.8 * (vdw_distance1 + vdw_distance2);
                        //if ( distance < vdwSum )
                        //    distance = vdwSum;
                        if ( distance < MIN_DISTANCE )
                            distance = MIN_DISTANCE;
                        double electrostatic = (charge1 * charge2 * COULOMB_CONSTANT) / (distance * distance);

                        // apply the combining rules for epsilon and sigma if necessary
                        double sigma = vdw_distance1;
                        double epsilon = vdw_depth1;

                        if ( atomClass1 != atomClass2 )
                            {
                                sigma = Math.sqrt(vdw_distance1 * vdw_distance2);
                                epsilon = Math.sqrt(vdw_depth1 * vdw_depth2);
                            }

                        // calculate the steric energy
                        double temp = Math.pow(sigma / distance, 6);
                        double steric = 4.0 * epsilon * temp * ( temp - 1.0 );

                        // return the result
                        energy += electrostatic + steric;
                    }
            }
        return energy;
    }

    public static String getDebugInteractionEnergy(Peptide peptide, Rotamer rotamer1, Rotamer rotamer2)
    {
        double energy = 0.0;

        TreeMap<Double,String> interactions = new TreeMap<>();
        for (Atom atom1 : rotamer1.atoms)
            {
                Integer atomType1 = atom1.type2;
                Integer atomClass1 = getOPLSClass(atomType1);

                double charge1 = getCharge(atomType1);
                double vdw_distance1 = getVDWDistance(atomClass1);
                double vdw_depth1 = getVDWDepth(atomClass1);
                
                for (Atom atom2 : rotamer2.atoms)
                    {
                        Integer atomType2 = atom2.type2;
                        Integer atomClass2 = getOPLSClass(atomType2);

                        double charge2 = getCharge(atomType2);
                        double vdw_distance2 = getVDWDistance(atomType2);
                        double vdw_depth2 = getVDWDepth(atomType2);

                        // calculate the coulombic energy
                        double distance = Vector3D.distance(atom1.position, atom2.position);
                        if ( distance < MIN_DISTANCE )
                            distance = MIN_DISTANCE;
                        double electrostatic = (charge1 * charge2 * COULOMB_CONSTANT) / distance;

                        // apply the combining rules for epsilon and sigma if necessary
                        double sigma = vdw_distance1;
                        double epsilon = vdw_depth1;

                        if ( atomClass1 != atomClass2 )
                            {
                                sigma = Math.sqrt(vdw_distance1 * vdw_distance2);
                                epsilon = Math.sqrt(vdw_depth1 * vdw_depth2);
                            }

                        // calculate the steric energy
                        double temp = Math.pow(sigma / distance, 6);
                        double steric = 4.0 * epsilon * temp * ( temp - 1.0 );

                        // return the result
                        energy += electrostatic + steric;

                        // make debug string
                        String debugString = String.format("%4s - %4s  dist = %5.2f   electrostatic: %12.2f  steric:  %12.2f  q1 = %7.4f  q2 = %7.4f  sigma = %7.4f  epsilon = %7.4f\n",
                                             peptide.getClosestAtomString(atom1), peptide.getClosestAtomString(atom2), distance,
                                             electrostatic, steric, charge1, charge2, epsilon, sigma);
                        interactions.put(electrostatic+steric, debugString);

                    }
            }

        String returnString = "";
        for (Double d : interactions.keySet())
            {
                returnString += interactions.get(d);
            }
        return returnString;
    }

    public static void main(String[] args)
    {
        DatabaseLoader.go();
        List<Peptide> sheets = BetaSheetGenerator.generateSheets(5, 5, 10000, 0.01);
        Peptide peptide = sheets.get(0);
        
        int sequenceLength = peptide.sequence.size();
        int forbiddenIndex = (sequenceLength/2) - 1;
        List<String> stringSequence = ImmutableList.of("standard_alanine", "arg", "asparagine", "aspartate", "glutamine", "glycine",
                                                       "histidine_hd", "isoleucine", "phenylalanine", "serine", "tryptophan", "standard_alanine");
        List<ProtoAminoAcid> protoAminoAcids = ProtoAminoAcidDatabase.getSpecificSequence(stringSequence);
        int j = 0;
        for (int i=0; i < sequenceLength; i++)
            {
                if ( i == forbiddenIndex || i == forbiddenIndex+1 )
                    continue;
                Residue residue = peptide.sequence.get(i);
                ProtoAminoAcid protoAminoAcid = protoAminoAcids.get(j);
                peptide = SidechainMutator.mutateSidechain(peptide, residue, protoAminoAcid);
                j++; 
            }
        new GaussianInputFile(peptide).write("test_peptides/test.gjf");
        List<Interaction> interactions = getInteractions(peptide);
        for (int i=0; i < 10; i++)
            System.out.println(interactions.get(i).toString(peptide));

        /*Atom atom1 = new Atom("H", new Vector3D(0.0, 0.0, 0.0), 28);
        Atom atom2 = new Atom("H", new Vector3D(0.0, 0.0, 2.0), 28);
        
        List<Atom> contents = ImmutableList.of(atom1, atom2);
        SimpleWeightedGraph<Atom,DefaultWeightedEdge> connectivity = new SimpleWeightedGraph<>(DefaultWeightedEdge.class);
        connectivity.addVertex(atom1);
        connectivity.addVertex(atom2);

        Molecule molecule = new Molecule("test molecule", contents, connectivity);
        System.out.println(molecule);
        
        List<Interaction> interactions = getInteractions(molecule);
        for (Interaction i : interactions)
            System.out.println(i.toString(molecule));*/
    }
}
