import java.io.*;
import java.util.*;
import java.util.concurrent.*;
import Jama.*;
import com.google.common.collect.*;
import org.apache.commons.math3.geometry.euclidean.threed.*;

/**
 * This class takes a bunch of peptides and aligns them all by C(alpha) position to minimize RMSD.
 * All structures are first translated to the origin.  Then a rigid body rotation is carried out using the
 * Kabsch algorithm to minimize the C(alpha) deviation between each structure and the anchor structure.
 */
public class Superposition implements Immutable, Serializable
{
    /** For serialization. */
    public static final long serialVersionUID = 1L;

    /** The coordinates of the C(alpha) atoms in the anchor structure. */
    public List<Vector3D> anchorGeometry;

    /**
     * The superimposed coordinates for each structure.  The outer index is the structure index.
     * The inner index is the residue number.  The residue numbering is assumed to be constant across 
     * all structures.  Only C(alpha) atoms are included.
     */
    public List<List<Vector3D>> resultGeometries;

    /** The root-mean-square Calpha deviations from the anchor structure. */
    public final List<Double> RMSDs;

    /** Standard constructor. */
    private Superposition(List<Vector3D> anchorGeometry, List<List<Vector3D>> resultGeometries, List<Double> RMSDs)
    {
        this.anchorGeometry = ImmutableList.copyOf(anchorGeometry);
        this.resultGeometries = ImmutableList.copyOf(resultGeometries);
        this.RMSDs = ImmutableList.copyOf(RMSDs);
    }

    /**
     * Static factory method for superimposing a list of peptides.
     * @param anchorPeptide all peptides will be aligned to this structure
     * @param peptides the peptides to align
     * @return the result of the superposition process
     */
    public static Superposition superimpose(Peptide anchorPeptide, List<Peptide> peptides)
    {
        int expectedNumberOfAtoms = anchorPeptide.contents.size();
        List<Vector3D> anchorPositions = center(extractCAs(anchorPeptide));
        List<List<Vector3D>> otherPositions = new ArrayList<>(peptides.size());
        int expectedSize = anchorPeptide.sequence.size();
        for (Peptide p : peptides)
            {
                if ( p.contents.size() != expectedNumberOfAtoms )
                    throw new IllegalArgumentException("unexpected number of atoms");
                if ( p.sequence.size() != expectedSize )
                    throw new IllegalArgumentException("sequence size mismatch");
                List<Vector3D> thesePositions = center(extractCAs(p));
                otherPositions.add(thesePositions);
            }
        List<List<Vector3D>> superimposedPositions = superimpose(anchorPositions, otherPositions);
        List<Double> RMSDs = new ArrayList<>(peptides.size());
        for (List<Vector3D> positions : superimposedPositions)
            RMSDs.add( getRMSD(anchorPositions,positions) );
        return new Superposition(anchorPositions, superimposedPositions, RMSDs);
    }

    /**
     * Static factory method for superimposing a list of peptides.
     * @param anchorPeptide all peptides will be aligned to this structure
     * @param structures a list of conformers of the anchorPeptide that will be aligned with the anchor peptide
     * @return the result of the superposition process
     */
    public static Superposition superimpose2(Peptide anchorPeptide, List<List<Vector3D>> structures)
    {
        int expectedNumberOfAtoms = anchorPeptide.contents.size();
        List<Vector3D> anchorPositions = center(extractCAs(anchorPeptide));
        List<Integer> anchorIndices = extractCAIndices(anchorPeptide);

        List<List<Vector3D>> otherPositions = new ArrayList<>(structures.size());
        int expectedSize = anchorPeptide.sequence.size();
        
        for (List<Vector3D> positions : structures)
            {
                if ( positions.size() != expectedNumberOfAtoms )
                    throw new IllegalArgumentException("unexpected number of atoms");
                List<Vector3D> currentCAPositions = new ArrayList<>(expectedSize);
                for (Integer i : anchorIndices)
                    currentCAPositions.add(positions.get(i));
                if ( currentCAPositions.size() != expectedSize )
                    throw new IllegalArgumentException("unexpected size");
                currentCAPositions = center(currentCAPositions);
                otherPositions.add(currentCAPositions);
            }
        List<List<Vector3D>> superimposedPositions = superimpose(anchorPositions, otherPositions);
        List<Double> RMSDs = new ArrayList<>(structures.size());
        for (List<Vector3D> positions : superimposedPositions)
            RMSDs.add( getRMSD(anchorPositions,positions) );
        return new Superposition(anchorPositions, superimposedPositions, RMSDs);
    }

    /**
     * Extracts the positions of all the alpha carbons in a peptide.
     * @param peptide the peptide to get the CA positions of
     * @return the CA positions in order of residue
     */
    public static List<Vector3D> extractCAs(Peptide peptide)
    {
        List<Vector3D> returnList = new ArrayList<>(peptide.sequence.size());
        for (Residue r : peptide.sequence)
            returnList.add( r.CA.position );
        return ImmutableList.copyOf(returnList);
    }

    /**
    * Extracts the atom indices of all the alpha carbons in a peptide.
    * @param peptide the peptide to get the CA indices of
    * @return the CA atom indices in order of residue
    */
    public static List<Integer> extractCAIndices(Peptide peptide)
    {
        List<Integer> returnList = new ArrayList<>(peptide.sequence.size());
        for (Residue r : peptide.sequence)
            returnList.add( peptide.contents.indexOf(r.CA ));
        return ImmutableList.copyOf(returnList);
    }

    /**
     * Centers the specified geometry at the origin.
     * @param positions the original positions
     * @return the centered positions
     */
    public static List<Vector3D> center(List<Vector3D> positions)
    {
        Vector3D centroid = Vector3D.ZERO;
        for (Vector3D v : positions)
            centroid = centroid.add(v);
        centroid = centroid.scalarMultiply(1.0/positions.size());
        List<Vector3D> newPositions = new ArrayList<>(positions.size());
        for (Vector3D v : positions)
            {
                Vector3D newVector = v.add( centroid.negate() );
                newPositions.add(newVector);
            }
        return ImmutableList.copyOf(newPositions);
    }

    /**
     * Performs a Kabsch superposition.  See <a href="http://en.wikipedia.org/wiki/Kabsch_algorithm">Kabsch algorithm</a>.
     * Positions are assumed to have been centered at the origin already.
     * @param anchorPositions the positions to align everything to
     * @param otherPositions lists of positions to align
     * @return the superimposed positions
     */
    public static List<List<Vector3D>> superimpose(List<Vector3D> anchorPositions, List<List<Vector3D>> otherPositions)
    {
        // form a matrix P that contains the coordinates of the anchor positions
        int size = anchorPositions.size();
        double[][] pre_matrixP = new double[size][3];
        for (int i=0; i < size; i++)
            {
                Vector3D v = anchorPositions.get(i);
                pre_matrixP[i][0] = v.getX();
                pre_matrixP[i][1] = v.getY();
                pre_matrixP[i][2] = v.getZ();
            }
        Matrix matrixP = new Matrix(pre_matrixP);
        Matrix matrixPtranspose = matrixP.transpose();

        // loop over all of the given sets of positions
        List<List<Vector3D>> returnList = new ArrayList<>(otherPositions.size());
        for (List<Vector3D> positions : otherPositions)
            {
                // form a matrix Q that contains the coordinates
                double[][] pre_matrixQ = new double[size][3];
                for (int i=0; i < size; i++)
                    {
                        Vector3D v = positions.get(i);
                        pre_matrixQ[i][0] = v.getX();
                        pre_matrixQ[i][1] = v.getY();
                        pre_matrixQ[i][2] = v.getZ();
                    }
                Matrix matrixQ = new Matrix(pre_matrixQ);

                // calculate the covariance matrix
                Matrix covarianceMatrix = matrixPtranspose.times(matrixQ);

                // get the singular valude decomposition
                SingularValueDecomposition SVD = new SingularValueDecomposition(covarianceMatrix);
                Matrix matrixV = SVD.getU();
                Matrix matrixW = SVD.getV();

                // check the sign of the rotation matrix
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

                // apply rotation
                List<Vector3D> newPositions = new ArrayList<>(positions.size());
                for (Vector3D v : positions)
                    {
                        currentVectorPreMatrix[0][0] = v.getX();
                        currentVectorPreMatrix[1][0] = v.getY();
                        currentVectorPreMatrix[2][0] = v.getZ();
                        Matrix currentVectorMatrix = new Matrix(currentVectorPreMatrix);
                        Matrix newVectorMatrix = rotationMatrix.times(currentVectorMatrix);
                        newVectorPreMatrix = newVectorMatrix.getArrayCopy();
                        Vector3D newVector = new Vector3D(newVectorPreMatrix[0][0], newVectorPreMatrix[1][0], newVectorPreMatrix[2][0]);
                        newPositions.add(newVector);
                    }

                // update the results
                returnList.add(ImmutableList.copyOf(newPositions));
            }
        return returnList;
    }

    /**
     * Computes the root-mean-square deviation between two lists of positions.
     * @param positions1 the first list of positions
     * @param positions2 the second list of positions
     * @return the RMSD
     */
    public static double getRMSD(List<Vector3D> positions1, List<Vector3D> positions2)
    {
        if (positions1.size() != positions2.size())
            throw new IllegalArgumentException("list size mismatch");
        double RMSD = 0.0; 
        for (int i=0; i < positions1.size(); i++)
            {
                Vector3D v1 = positions1.get(i);
                Vector3D v2 = positions2.get(i);
                double distance = Vector3D.distance(v1,v2);
                RMSD += distance * distance;
            }
        RMSD = Math.sqrt(RMSD / positions1.size());
        return RMSD;
    }

    @Override
    public String toString()
    {
        String returnString = String.format("%d superimposed geometries\n", RMSDs.size());
        for (Double d : RMSDs)
            returnString += String.format("%.4f  ", d);
        return returnString;
    }

    @Override
    public boolean equals(Object obj)
    {
      if ( obj == null )
            return false;
        if ( obj == this )
            return true;
        if ( !(obj instanceof Superposition) )
            return false;

        Superposition s = (Superposition)obj;
        if ( s.RMSDs.equals(RMSDs) )
            return true;
        return false;
    }

    @Override
    public int hashCode()
    {
        return Objects.hash(RMSDs);
    }
}
