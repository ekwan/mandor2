import java.util.*;
import java.io.*;
import com.google.common.collect.*;
import org.apache.commons.math3.geometry.euclidean.threed.*;

/**
 * This class calculates the solvent-exposed surface area of a molecule using numerical methods.
 * It produces around each atom a mesh of evenly-distributed points, and then it checks for each 
 * point whether or not it's contained in the sphere around any other atomic center.  If the point
 * is exposed to solvent, then it is added to a tally.  The surface area corresponding to each
 * point is calculated, and therefore we can determine the total surface area accessible to each 
 * atom.
 *
 * There is an optional probe radius (1.40 is typical, as it is the radius of water).  This means
 * that we get the area accessible to rolling a ball of that radius around the molecule.  This feature
 * is untested.
 *
 * The grid is generated using the golden section method.  The points are generated in cylindrical
 * coordinates (angle, radius, height).  If we want n points, let i run from 0 to n-1 (such that
 * we'll have a point at zero angle):
 *
 * angle = GA * i
 * radius = sqrt(1 - z_i^2)
 * z_i = (1 - 1/n) * (1 - 2i/n-1)
 *
 * where GA is the golden angle: pi - 3*sqrt(5).  It's pretty amazing that works!
 * 
 * Mostly Greg and I just followed the recipe here:
 *
 * http://blog.marmakoide.org/?p=1
 * 
 * The Shrake-Rupley algorithm is described everywhere, but this is a good link:
 *
 * http://boscoh.com/protein/calculating-the-solvent-accessible-surface-area-asa.html
 *
 * It all makes sense except for the probe radius part, which we find suspicious.
 */
public class ShrakeRupleyCalculator extends SurfaceAreaCalculator implements Immutable
{
    /** The database of atomic radii in angstroms. */
    public static final Map<Element,Double> RADII;

    /** The golden angle.  Very badly approximable by rationals.  */
    public static final double GA = Math.PI*(3.0-Math.sqrt(5.0));

    /** Static initializer. */
    static
        {
            Map<Element,Double> tempMap = new HashMap<>();
            
            // these are the BONDII radii
            tempMap.put(Element.CARBON, 1.70);
            tempMap.put(Element.NITROGEN, 1.55);
            tempMap.put(Element.OXYGEN, 1.52);
            tempMap.put(Element.HYDROGEN, 1.20);
            tempMap.put(Element.SULFUR, 1.80);

            // these are the LCPO radii
            /*tempMap.put(Element.CARBON, 1.70);
            tempMap.put(Element.NITROGEN, 1.65);
            tempMap.put(Element.OXYGEN, 1.60);
            tempMap.put(Element.HYDROGEN, 0.00);  // treated implicitly
            tempMap.put(Element.SULFUR, 1.90);*/

            RADII = ImmutableMap.copyOf(tempMap);
        }

    /** Number of points to use in the mesh. */
    public final int meshSize;

    /** The points in the mesh. */
    private final Vector3D[] grid;

    /** Probe radius to use. */
    public final double probeRadius;

    public static final ShrakeRupleyCalculator INSTANCE = new ShrakeRupleyCalculator(1000, 0.0);

    /**
     * Creates a surface area calculator.
     */
    public ShrakeRupleyCalculator(int meshSize, double probeRadius)
    {
        this.meshSize = meshSize;
        this.probeRadius = probeRadius;
        this.grid = pointMesh();
    }

    /**
     * Calculates the solvent accessible surface area of the specified molecule.
     * @param molecule the molecule to analyze
     * @return the SASA by atom in angstroms^2
     */
    @Override
    public List<Double> calculateSASA(Molecule molecule)
    {
        // get the centers for every atom
        int numberOfAtoms = molecule.contents.size();
        List<Vector3D> centers = new ArrayList<>(numberOfAtoms);
        for (Atom a : molecule.contents)
            centers.add(a.position);

        // get radii for every atom
        List<Double> radii = getRadii(molecule);

        // get mesh of points for every atom
        List<Vector3D[]> meshes = getMeshes(centers, radii);
        //System.out.println(Arrays.toString(meshes.get(0)));

        // for every atom, consider every point in its mesh
        List<Double> SASAlist = new ArrayList<>(numberOfAtoms);
        for (int i=0; i < numberOfAtoms; i++)
            {
                Vector3D center1 = centers.get(i);
                double radius1 = radii.get(i);
                Vector3D[] thisMesh = meshes.get(i);
                int exposed = 0;
                for (Vector3D meshPoint : thisMesh)
                    {
                        // consider all atoms that are adjacent to this one
                        boolean buried = false;
                        for (int j=0; j < numberOfAtoms; j++)
                            {
                                // skip if this is the same atom
                                if ( i==j )
                                    continue;
                                
                                // ignore atoms that are too far away to interact
                                Vector3D center2 = centers.get(j);
                                double radius2 = radii.get(j);
                                double distance = Vector3D.distance(center1,center2);
                                if ( distance > radius1 + radius2 )
                                    continue;

                                // check if this atom is buried
                                if ( isBuried(meshPoint, center2, radius2, probeRadius) )
                                    {
                                        buried = true;
                                        break;
                                    }
                            }
                        if ( !buried )
                            exposed++;
                    }
               
                // compute the surface area per point for this atom
                double totalSurfaceArea = 4.0 * Math.PI * radius1 * radius1;
                double areaPerPoint = totalSurfaceArea / meshSize;

                // compute the surface area for this atom
                double SASA = areaPerPoint * exposed;
                SASAlist.add(SASA);
            }
        
        // return the result
        return ImmutableList.copyOf(SASAlist);
    }

    /**
     * Produce a list of points in three dimensions lying in 
     * the sphere centered at the origin.  
     * @param N The number of points in the list. 
     * @return A list of N points lying on the sphere centered at the origin.  
     */
    private Vector3D[] pointMesh(int N) {
        Vector3D[] output = new Vector3D[N];
        double angle = 0.0d;
        double z = (1-1/((double)N));
        double r = Math.sqrt(1-z*z);
        double deltaZ = 2*z/((double)N-1.0d);
        for (int i = 0; i < N; i++) {
            output[i] = new Vector3D(r*Math.cos(angle),r*Math.sin(angle),z);
            angle += GA;
            z -= deltaZ;
            r = Math.sqrt(1-z*z);
        }
        return output;
    }

    /**
     * Same as above, but set the number of points to a constant value.  
     */
    private Vector3D[] pointMesh() {
        return pointMesh(meshSize);
    }

    /**
     * Transform a list of points in three dimensions by multiplying them 
     * them by a common scalar and then shifting them all by a common 
     * translation vector.  Do not change the underlying list; just create 
     * a new one.  
     * @param s The scalar multiple.  
     * @param v The translation vector.  
     * @return The list of transformed points.  
     */
    private Vector3D[] transform(Vector3D[] input, double s, Vector3D v) {
        Vector3D[] output = new Vector3D[input.length];
        for (int i = 0; i < output.length; i++) {
            output[i] = input[i].scalarMultiply(s).add(v);
        }
        return output;
    }

    /**
     * Produce a list of arrays of points.  Each array contains point lying on 
     * a common sphere with a given center and radius.  
     * @param centers The centers of the spheres.  
     * @param radii The radii of the spheres.  
     * @return The list of arrays of points.  
     */
    private List<Vector3D[]> getMeshes(List<Vector3D> centers, List<Double> radii) {
        List<Vector3D[]> output = new LinkedList<>();
        Vector3D[] baseMesh = grid;
        for (int i = 0; i < centers.size(); i++) {
            output.add(transform(baseMesh, radii.get(i) + probeRadius, centers.get(i)));
        }
        return output;
    }

    /**
     * Returns the atomic radii for the given molecule.
     * @param molecule the molecule to analyze
     * @return the radii in angstroms ordered by atom index
     */
    public List<Double> getRadii(Molecule molecule)
    {
        List<Double> returnList = new ArrayList<>(molecule.contents.size());
        for (Atom a : molecule.contents)
            returnList.add(RADII.get(a.element));
        return returnList;
    }

    /**
     * Tests whether the specified point is buried in the specified atom.
     * @param point the candidate point
     * @param atomCenter the center of the atom in question
     * @param atomRadius the radius in angstroms of the atom in question
     * @param probeRadius the radius of the probe
     * @return true if the point is buried
     */
    public static boolean isBuried(Vector3D point, Vector3D atomCenter, double atomRadius, double probeRadius)
    {
        double distance = Vector3D.distance(point, atomCenter);
        if ( distance < atomRadius + probeRadius )
            return true;
        return false;
    }

    /** For testing. */
    public static void main(String[] args)
    {
        TinkerXYZOutputFile testFile = new TinkerXYZOutputFile("amino_acids/hairpin_minimized.xyz");
        ShrakeRupleyCalculator calculator = new ShrakeRupleyCalculator(100, 0.00);
        List<Double> SASAlist = null;
        long startTime = System.currentTimeMillis();
        for (int i=0; i < 100; i++)
            {
                SASAlist = calculator.calculateSASA(testFile.molecule);
                //System.out.printf("%d \r", i);
            }
        long endTime = System.currentTimeMillis();
        long difference = endTime - startTime;
        System.out.println("\n" + difference/100);
        
        double sum = 0.0;
        for (int i=0; i < SASAlist.size(); i++)
            {
                Double d = SASAlist.get(i);
                sum += d;
                //System.out.printf("%3d  %.2f \n", i+1, d);
            }
        System.out.printf("\ntotal SASA = %.4f angstroms^2\n", sum);
        
    }
} // end of class ShrakeRupleyCalculator 
