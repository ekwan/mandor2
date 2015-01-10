import java.io.*;
import java.util.*;
import org.apache.commons.math3.geometry.euclidean.threed.*;
import org.jgrapht.*;
import org.jgrapht.graph.*;
import org.jgrapht.alg.*;
import com.google.common.collect.*;

/** A class representing Tinker XYZ Files in disk. Allows for reading in Molecules from xyz files **/
public class TinkerXYZOutputFile extends OutputFileFormat implements Immutable {

    /** The Molecule that the XYZFile describes. */
    public final Molecule molecule;
    
    /**
     * Reads TINKER geometry, connectivity, and atom types.
     * @param fileName the location of the XYZ file
    */
    public TinkerXYZOutputFile(String fileName)
    {
	    super(fileName);
	
        // read in atoms, name, and connectivity
	    String name = "";
	    List<Atom> atoms = new ArrayList<>();
	    SimpleWeightedGraph<Atom,DefaultWeightedEdge> connectivity = new SimpleWeightedGraph<>(DefaultWeightedEdge.class);

	    // holds connections till end of loop so that all atoms needed for edges have already been added to atoms
	    List<List<Integer>> connections = new ArrayList<>();
	
	    boolean isFirstLine = true;
        for (List<String> line : fileContents)
            {
                if (isFirstLine)
                    {
                        name = line.get(1);
                        isFirstLine = false;
                    }
                else
                    {
                        Vector3D position = new Vector3D(Double.parseDouble(line.get(2)),
                                                         Double.parseDouble(line.get(3)),
                                                         Double.parseDouble(line.get(4)));
                        int type1 = Integer.parseInt(line.get(5));

                        // only sets AMOEBA atom type
                        // OPLS atom types and surface tension are set to zero
			            Atom newAtom = new Atom(line.get(1), position, type1, 0, 0.0);
			            atoms.add(newAtom);
		
			            List<Integer> currentConnections = new ArrayList<>();
			            for (int i = 6; i<line.size() ;i++)
				            currentConnections.add(Integer.parseInt(line.get(i)));

			            connections.add(currentConnections);
			            connectivity.addVertex(newAtom);
		            }
	        }
	
        // build connectivity graph
	    int currentAtomIndex = 1;
	    for (List<Integer> toAtomIndices : connections)
	        {
		        // the atom from which we will be making connections ; the index in connections -1 is this atom's number in the XYZ File
		        Atom fromAtom = atoms.get(currentAtomIndex -1 );
		
		        // read in remaining atoms with bond to atom1 and their bond orders
		        for (Integer toAtomIndex : toAtomIndices)
		            {
			            Atom toAtom = atoms.get(toAtomIndex - 1);
			            
                        // check that bond we are trying to create does not already exist
			            if (connectivity.getEdge(fromAtom, toAtom) == null)
			                {
				                // create an edge in the connectivity graph
                                DefaultWeightedEdge thisEdge = connectivity.addEdge(fromAtom, toAtom);
				
                                // TINKER does not recognize bond orders, so set all to 1
                                double bondOrder = 1.0;
                                connectivity.setEdgeWeight(thisEdge, bondOrder);
                            }
                    }
		        currentAtomIndex++;
	    }

	    molecule = new Molecule(name, atoms, connectivity);
    }
    
    public static void main(String args[])
    {
        TinkerXYZOutputFile file = new TinkerXYZOutputFile("amino_acids/Ala.xyz");
        System.out.println(file.molecule);
    }

}
