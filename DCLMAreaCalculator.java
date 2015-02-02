import java.util.*;
import java.io.*;
import com.google.common.collect.*;
import org.apache.commons.math3.geometry.euclidean.threed.*;
import java.util.concurrent.*;

/**
 * This class calculates the solvent-exposed surface area of a molecule using numerical methods.
 * This uses the DCLM method reported in J. Comp. Chem. 1995, 16(3), 273-284.  This is a more
 * efficient version of Shrake-Rupley that does not require any special parameters beyond the
 * radii of thea toms themselves.  The process of neighbor list generation is improved by only
 * checking atom neighbors on a cubic grid where the spacing is 2*r_max, where r_max is the
 * largest radius.  That means that when two atoms are in contact, so too must their grid cells.
 *
 * The area calculation is improved as well.  In the Shrake-Rupley, we check if all points p
 * on the surface of atom i are enclosed in neighboring atoms j.  Points are checked regardless
 * of spatial proximity.  Now, we first subdivide the volume occupied by i and j into cubic
 * boxes.  We look for overlaps in the boxes, and then only consider points i that are inside
 * the boxes.
 *
 * One instance per conformation should be used; this class is not thread safe.
 */
public class DCLMAreaCalculator extends SurfaceAreaCalculator implements Immutable
{
    /** The database of atomic radii in angstroms. */
    public static final Map<Element,Double> RADII;

    /** The maximum of the atomic radii in angstroms. */
    public static final double MAX_RADIUS;

    /** The golden angle.  Very badly approximable by rationals.  */
    public static final double GA = Math.PI*(3.0-Math.sqrt(5.0));

    /** A 3d lattice of lists of atom centers and radii.  */
    private Map<Vector3D,Double>[][][] atomGrid;

    /** Number of points to use in the mesh. */
    public static final int MESH_SIZE = 1000;

    /** Number of subdivisions of the unit cube along each axis. Used to partition the mesh points.  */
    public static final int M;

    /** A non-final version of M.  We set the value of this in makeGrid and then set M = mmm.  */
    public static int mmm = 0;

    /** Points distributed approximately evenly over the surface of the unit sphere.  */
    private static final Vector3D[] MESH;

    /** The same points as in MESH, but organized into boxes indexed by their bottom corners.  */
    private static final Map<Vector3D,List<Vector3D>> GRID;

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

            RADII = ImmutableMap.copyOf(tempMap);

            // find the maximum of the atomic radii
            double max = -1.0d;
            for (Element e : RADII.keySet())
                if (RADII.get(e) > max) max = RADII.get(e);
            MAX_RADIUS = max;

            MESH = getMesh(MESH_SIZE);
            GRID = makeGrid(MESH);
            M = mmm;
        }

    /** Probe radius to use. */
    public final double probeRadius;

    /**
     * Creates a surface area calculator.  The probe radius is currently not supported and is ignored.
     */
    public DCLMAreaCalculator(double probeRadius)
    {
        this.probeRadius = probeRadius;
        this.atomGrid = null;
    }

    /**
     * Calculates the solvent accessible surface area of the specified molecule.
     * @param molecule the molecule to analyze
     * @return the SASA by atom in angstroms^2
     */
    @Override
    @SuppressWarnings("unchecked")
    public List<Double> calculateSASA(Molecule molecule)
    {
        double width = 2*MAX_RADIUS;
        // get the centers for every atom
        // while we're at it, determine min and max values for the x-, 
        // y-, and z-coordinates of atom centers.  
        double xmax = 0.0d;
        double xmin = 0.0d;
        double ymax = 0.0d;
        double ymin = 0.0d;
        double zmax = 0.0d;
        double zmin = 0.0d;
        Vector3D pos = new Vector3D(0.0d,0.0d,0.0d);
        boolean first = true;
        double x = 0.0d;
        double y = 0.0d;
        double z = 0.0d;
        int numberOfAtoms = molecule.contents.size();
        List<Vector3D> centers = new ArrayList<>(numberOfAtoms);
        for (Atom a : molecule.contents) {
            pos = a.position;
            x = pos.getX();
            y = pos.getY();
            z = pos.getZ();
            if (first) {
                xmax = x;
                xmin = x;
                ymax = y;
                ymin = y;
                zmax = z;
                zmin = z;
                first = false;
            } else {
                if (x > xmax) xmax = x;
                if (x < xmin) xmin = x;
                if (y > ymax) ymax = y;
                if (y < ymin) ymin = y;
                if (z > zmax) zmax = z;
                if (z < zmin) zmin = z;
            }
            centers.add(pos);
        }

        // get radii for every atom
        List<Double> radii = getRadii(molecule);

        // find the number of boxes in the atomGrid
        int numX = (int)((xmax-xmin)/width) + 1;
        int numY = (int)((ymax-ymin)/width) + 1;
        int numZ = (int)((zmax-zmin)/width) + 1;

        // fill the atomGrid
        atomGrid = (Map<Vector3D,Double>[][][])new Map<?,?>[numX][numY][numZ];
        for (int i=0; i < numX; i++)
            for (int j=0; j < numY; j++)
                for (int k=0; k < numZ; k++)
                    atomGrid[i][j][k] = new HashMap<Vector3D,Double>();
        for (int i = 0; i < centers.size(); i++) {
            pos = centers.get(i);
            int xIndex = (int)((pos.getX()-xmin)/width);
            int yIndex = (int)((pos.getY()-ymin)/width);
            int zIndex = (int)((pos.getZ()-zmin)/width);
            atomGrid[xIndex][yIndex][zIndex].put(pos, radii.get(i));
        }

        // calculate surface area
        List<Double> SASAlist = new ArrayList<>(numberOfAtoms);
        for (int l = 0; l < numberOfAtoms; l++) SASAlist.add(0.0d);
        // map from centers of neighbors of a given central atom to their radii
        Map<Vector3D,Double> nbhd = null;
        // current atomic radius
        double r = 0.0d;
        // the number of occluded mesh points on the current central atom
        int occluded = 0;
        // iterate over all grid regions, and all atoms
        for (int i = 0; i < numX; i++) {
            for (int j = 0; j < numY; j++) {
                for (int k = 0; k < numZ; k++) {
                    // skip empty sectors
                    if (atomGrid[i][j][k].isEmpty()) continue;
                    for (Vector3D v : atomGrid[i][j][k].keySet()) {
                        r = atomGrid[i][j][k].get(v);
                        nbhd = getNeighbors(i,j,k,v,r);
                        nbhd = transform(nbhd,v,r);
                        occluded = occluded(nbhd);
                        if ( occluded > 0 )
                            SASAlist.set(centers.indexOf(v),4.0d*Math.PI*r*r*(MESH_SIZE-occluded)/MESH_SIZE);
                    }
                }
            }
        }
        // return the result
        return ImmutableList.copyOf(SASAlist);
    } // end of method calculateSASA

    /**
     * Produce a list of points in three dimensions lying in 
     * the sphere centered at the origin.  
     * @param N The number of points in the list. 
     * @return A list of N points lying on the sphere centered at the origin.  
     */
    public static Vector3D[] getMesh(int N) {
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
     * Takes a set of points lying on the surface of the sphere.  A cube is inscribed around the sphere
     * and chopped up into n^3 sub-cubes.  Points on the surface of the sphere will be assigned to the
     * sub-cubes.
     * The keys represent the lowest xyz corner of each sub-cube.
     * The values represent the surface points lying in each sub-cube.
     */
    public static Map<Vector3D,List<Vector3D>> makeGrid(Vector3D[] grid)
    {
        int n = grid.length;
        int cubeSize = (int)Math.floor(Math.cbrt(n/2));
        // a kludge so that we can assign a static final variable using quantities from this method
        mmm = cubeSize;
        double cubeSpacing = 2.0 / cubeSize;
        List<Vector3D> keys = new ArrayList<>((int)Math.pow(cubeSize,3));
        for (double i=-1.0; i < 1.0-cubeSpacing; i+=cubeSpacing)
            {
                for (double j=-1.0; j < 1.0-cubeSpacing; j+=cubeSpacing)
                    {
                        for (double k=-1.0; k < 1.0-cubeSpacing; k+=cubeSpacing)
                            {
                                Vector3D v = new Vector3D(i,j,k);
                                keys.add(v);
                            }
                    }
            }
        List<Vector3D> assigned = new ArrayList<>();
        Map<Vector3D,List<Vector3D>> tempMap = new HashMap<>();
        for (Vector3D minVector : keys)
            {
                // determine the bounds of this cube
                double minX = minVector.getX();
                double minY = minVector.getY();
                double minZ = minVector.getZ();
                double maxX = minX + cubeSpacing;
                double maxY = minY + cubeSpacing;
                double maxZ = minZ + cubeSpacing;
                List<Vector3D> tempList = new ArrayList<>();
                for (Vector3D v : grid)
                    {
                        double x = v.getX();
                        double y = v.getY();
                        double z = v.getZ();
                        if ( minX <= x && x <= maxX &&
                             minY <= y && y <= maxY &&
                             minZ <= z && z <= maxZ    )
                            {
                                // check for double assignments
                                if ( assigned.contains(v) )
                                    throw new IllegalArgumentException("already assigned");

                                // this vector is within this sub-cube
                                assigned.add(v);
                                tempList.add(v);
                            }
                    }
                if (tempList.size() > 0)
                    {
                        tempMap.put(minVector, ImmutableList.copyOf(tempList));
                    }
            }
        if ( assigned.size() != grid.length )
            throw new IllegalArgumentException("not everything got assigned");
        return ImmutableMap.copyOf(tempMap);
    }

    /**
     * Locates the intersection of a sphere with the unit sphere and returns the minimum and
     * maximum coordinates of intersection along an axis.  This is x_l and x_u in Figure 3 of
     * the paper.
     * @param center the center of the sphere
     * @param radius the radius of the sphere
     * @return [x_l, x_u, y_l, y_u, z_l, z_u]
     */
    public static double[] getProjection(Vector3D center, double radius)
    {
        // check invariants
        if ( center == null || radius <= 0.0 )
            throw new IllegalArgumentException("unexpected center or radius");
    
        // calculate some quantities we need
        double d_ij_squared = center.getNormSq();
        double d_ij = center.getNorm(); // we need an extra square root here not mentioned in the paper
        double r_squared = radius * radius;
        double C = (d_ij_squared + 1.0 - r_squared) / (2.0d * d_ij_squared);
        double D = (1.0 / d_ij_squared) - C*C;
        double x_j = center.getX();
        double x_j_squared = x_j * x_j;
        double y_j = center.getY();
        double y_j_squared = y_j * y_j;
        double z_j = center.getZ();
        double z_j_squared = z_j * z_j;

        // make the return array
        double[] returnArray = new double[6];

        // compute the projection onto the x-axis
        double cos_cos = x_j * C;
        double sin_sin = Math.sqrt(D * (y_j_squared + z_j_squared));
        double cos_alpha = x_j / d_ij;
        double cos_omega = (cos_alpha != 0.0d) ? (cos_cos / cos_alpha) : 0.0d;
        if (cos_alpha == 0.0d) {
            if (d_ij_squared + 1 <= r_squared) {
                returnArray[0] = -1.0d;
                returnArray[1] =  1.0d;
            } else {
                returnArray[0] = cos_cos - sin_sin;
                returnArray[1] = cos_cos + sin_sin;
            }
        } else {
            returnArray[0] = (cos_alpha <= -cos_omega) ? -1.0d : cos_cos - sin_sin;
            returnArray[1] = (cos_alpha >=  cos_omega) ?  1.0d : cos_cos + sin_sin;
        }
        // compute the projection onto the y-axis
        cos_cos = y_j * C;
        sin_sin = Math.sqrt(D * (x_j_squared + z_j_squared));
        cos_alpha = y_j / d_ij;
        cos_omega = (cos_alpha != 0.0d) ? (cos_cos / cos_alpha) : 0.0d;
        if (cos_alpha == 0.0d) {
            if (d_ij_squared + 1 <= r_squared) {
                returnArray[2] = -1.0d;
                returnArray[3] =  1.0d;
            } else {
                returnArray[2] = cos_cos - sin_sin;
                returnArray[3] = cos_cos + sin_sin;
            }
        } else {
            returnArray[2] = (cos_alpha <= -cos_omega) ? -1.0d : cos_cos - sin_sin;
            returnArray[3] = (cos_alpha >=  cos_omega) ?  1.0d : cos_cos + sin_sin;
        }
        // compute the projection onto the z-axis
        cos_cos = z_j * C;
        sin_sin = Math.sqrt(D * (y_j_squared + x_j_squared));
        cos_alpha = z_j / d_ij;
        cos_omega = (cos_alpha != 0.0d) ? (cos_cos / cos_alpha) : 0.0d;
        if (cos_alpha == 0.0d) {
            if (d_ij_squared + 1 <= r_squared) {
                returnArray[4] = -1.0d;
                returnArray[5] =  1.0d;
            } else {
                returnArray[4] = cos_cos - sin_sin;
                returnArray[5] = cos_cos + sin_sin;
            }
        } else {
            returnArray[4] = (cos_alpha <= -cos_omega) ? -1.0d : cos_cos - sin_sin;
            returnArray[5] = (cos_alpha >=  cos_omega) ?  1.0d : cos_cos + sin_sin;
        }

        // all done!
        return returnArray;
    } // end of method getProjection

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
     * Scale and translate a collection of atoms.  
     * The intention is to put the neighbors of a given central atom in
     * relative coordinates.  
     * @param neighbors A map from centers of neighboring atoms to their radii.  
     * @param center The center of the given central atom.  
     * @param radius The radius of the given central atom.  
     * @return The same neighboring atoms expressed in relative coordinates.  
     */
    public Map<Vector3D,Double> transform(Map<Vector3D,Double> neighbors,Vector3D center,double radius)
    {
        double div = 1.0d/radius;
        Map<Vector3D,Double> output = new HashMap<Vector3D,Double>();
        for (Vector3D c : neighbors.keySet()) {
            output.put(c.subtract(center).scalarMultiply(div),neighbors.get(c)*div);
        }
        return output;
    }

    /**
     * Returns the neighbors of a given atom.  
     * @param center The position of the atom whose neighbors we wish to find.  
     * @param x The x-coordinate of the grid sector in which the points lies.  
     * @param y The y-coordinate of the grid sector in which the points lies.  
     * @param z The z-coordinate of the grid sector in which the points lies.  
     * @param rad The radius of the atom whose neighbors we wish to find.  
     * @return A Map from positions of centers of neighboring atoms to their radii.  
     */
    private Map<Vector3D,Double> getNeighbors(int x,int y,int z,Vector3D center,double rad) {
        Map<Vector3D,Double> output = new HashMap<Vector3D,Double>();
        for (int i = x-1; i < x+2; i++) {
            if (i<0||i>=atomGrid.length) continue;
            for (int j = y-1; j < y+2; j++) {
                if (j<0||j>=atomGrid[i].length) continue;
                for (int k = z-1; k < z+2; k++) {
                    if (k<0||k>=atomGrid[i][j].length) continue;
                    Map<Vector3D,Double> current = atomGrid[i][j][k];
                    for (Vector3D c : current.keySet()) {
                        if (c.equals(center)) continue;
                        if (c.distance(center)<rad+current.get(c)) output.put(c,current.get(c));
                    }
                }
            }
        }
        return output;
    }

    /**
     * Returns the number of mesh points occluded by a given set of neighbors.  
     * @param neighbors A map of atom centers and radii for the neighboring atoms.  
     * @return The number of mesh points occluded by the given neighboring atoms.  
     */
    private int occluded(Map<Vector3D,Double> neighbors) {
        // for each sector, make a list of neighboring atoms that might overlap it
        Map<Vector3D,List<Vector3D>> overlaps = new HashMap<>();
        for (Vector3D v : GRID.keySet()) {
            overlaps.put(v,new ArrayList<Vector3D>());
        }
        double[] bounds = null;
        double del = 2.0d / M;
        for (Vector3D v : neighbors.keySet()) {
            bounds = getProjection(v,neighbors.get(v));
            for (Vector3D w : overlaps.keySet()) {
                // on each axis check that the interval between the coord of w
                // and the coord of w plus del does not end before or begin
                // after the interval spanned by v
                if (!(w.getX() > bounds[1] || w.getX()+del < bounds[0]) &&
                    !(w.getY() > bounds[3] || w.getY()+del < bounds[2]) &&
                    !(w.getZ() > bounds[5] || w.getZ()+del < bounds[4]))
                {
                    overlaps.get(w).add(v);
                }
            }
        }

        // now iterate through all points and check to see if they're occluded
        int output = 0;
        // if there is a neighbor that occludes something, remember it
        Vector3D match = null;
        double matchRadSq = 0.0d;
        for (Vector3D v : GRID.keySet()) {
            // don't check these points if there are no neighbors overlapping this sector
            if (overlaps.get(v).isEmpty()) continue;
            // iterate over all mesh points in the sector
            meshloop:
            for (Vector3D m : GRID.get(v)) {
                if (match != null && match.distanceSq(m) < matchRadSq) {
                    output++;
                    continue;
                }
                for (Vector3D nbr : overlaps.get(v)) {
                    if (!nbr.equals(match) && nbr.distanceSq(m) < (matchRadSq = neighbors.get(nbr)*neighbors.get(nbr))) {
                        match = nbr;
                        output++;
                        continue meshloop;
                    }
                }
                match = null;
            }
        }
        return output;
    }

    /** For testing. */
    public static void main(String[] args)
    {
        TinkerXYZOutputFile testFile = new TinkerXYZOutputFile("amino_acids/hairpin_minimized.xyz");
        ShrakeRupleyCalculator calculator1 = new ShrakeRupleyCalculator(10000,0.00);
        DCLMAreaCalculator calculator2 = new DCLMAreaCalculator(0.00);
        System.out.print("Calculating Shrake-Rupley area...");
        List<Double> SASAlist1 = calculator1.calculateSASA(testFile.molecule);
        System.out.println("done.");
        List<Double> SASAlist2 = null;
        long startTime = System.currentTimeMillis();
        for (int i=0; i < 100; i++)
            {
                SASAlist2 = new DCLMAreaCalculator(0.0).calculateSASA(testFile.molecule);
                //System.out.printf("%d \r", i);
            }
        long endTime = System.currentTimeMillis();
        long difference = endTime - startTime;
        System.out.println("\n" + difference/100);        
        double sum1 = 0.0;
        double sum2 = 0.0;
        int x1 = 0;
        int x2 = 0;
        for (int i=0; i < testFile.molecule.contents.size(); i++)
            {
                sum1 += SASAlist1.get(i);
                sum2 += SASAlist2.get(i);
                double diff = (SASAlist2.get(i) - SASAlist1.get(i))*100.0/SASAlist1.get(i);
                System.out.printf("%4d %7.4f %7.4f %5.1f%%\n", i+1, SASAlist1.get(i), SASAlist2.get(i), diff);
            }
        System.out.printf("Shrake-Rupley area (A^2): %.4f\n", sum1);
        System.out.printf("DCLM area (A^2):          %.4f\n", sum2);
    }

} // end of class DCLMAreaCalculator 
