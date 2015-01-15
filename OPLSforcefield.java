import java.io.*;
import java.util.*;
import com.google.common.collect.*;

/**
 * This class reads parameters from the OPLS force field file and stores them in various maps.
 */
public class OPLSforcefield implements Singleton
{
    /** Map from OPLS atom types to atom classes. */
    public static final Map<Integer,Integer> CLASS_MAP;

    /** Map from atom classes to VDW distances (sigma) in A. */
    public static final Map<Integer,Double> VDW_DISTANCE_MAP;

    /** Map from atom classes to VDW depths (epsilon) in kcal/mol. */
    public static final Map<Integer,Double> VDW_DEPTH_MAP;

    /** Map from atom types to partial charges. */
    public static final Map<Integer,Double> CHARGE_MAP;

    private OPLSforcefield()
    {
        throw new IllegalArgumentException("not instantiable");
    }

    static
        {
            // make temporary fields
            Map<Integer,Integer> tempClassMap = new HashMap<>();
            Map<Integer,Double> tempDistanceMap = new HashMap<>();
            Map<Integer,Double> tempDepthMap = new HashMap<>();
            Map<Integer,Double> tempChargeMap = new HashMap<>();

            // read forcefield file
            OutputFileFormat forcefieldFile = new OutputFileFormat(Settings.OPLS_FORCEFIELD_FILENAME) {};

            for (List<String> fields : forcefieldFile.fileContents)
                {
                    if ( fields.size() == 0 )
                        continue;
                    if      ( fields.get(0).equals("atom")  && fields.size() >= 3 )
                        {
                            Integer typeNumber = Integer.valueOf(fields.get(1));
                            Integer classNumber = Integer.valueOf(fields.get(2));
                            if ( tempClassMap.containsKey(typeNumber) )
                                throw new IllegalArgumentException("duplicate type/class line\n" + fields.toString());
                            tempClassMap.put(typeNumber, classNumber);
                        }
                    else if ( fields.get(0).equals("vdw")   && fields.size() >= 4 )
                        {
                            Integer classNumber = Integer.valueOf(fields.get(1));
                            Double distance     = Double.valueOf(fields.get(2));
                            Double depth        = Double.valueOf(fields.get(3));
                            if ( tempDistanceMap.containsKey(classNumber) || tempDepthMap.containsKey(classNumber) )
                                throw new IllegalArgumentException("duplicate vdw line\n" + fields.toString());
                            tempDistanceMap.put(classNumber, distance);
                            tempDepthMap.put(classNumber, depth);
                        }
                    else if ( fields.get(0).equals("charge") && fields.size() >= 3 )
                        {
                            Integer type   = Integer.valueOf(fields.get(1));
                            Double  charge = Double.valueOf(fields.get(2));
                            if ( tempChargeMap.containsKey(type) )
                                throw new IllegalArgumentException("duplicate charge line\n" + fields.toString());
                            tempChargeMap.put(type,charge);
                        }
                }

            // set permanent fields
            CLASS_MAP = ImmutableMap.copyOf(tempClassMap);
            VDW_DISTANCE_MAP = ImmutableMap.copyOf(tempDistanceMap);
            VDW_DEPTH_MAP = ImmutableMap.copyOf(tempDepthMap);
            CHARGE_MAP = ImmutableMap.copyOf(tempChargeMap);
        }

    /** Forces the database to load. */
    public static void load()
    {
        System.out.printf("%d atom types have been loaded.\n", CLASS_MAP.size());
        System.out.printf("%d vdw classes have been loaded.\n", VDW_DISTANCE_MAP.size());
        System.out.printf("%d charge parameters have been loaded.\n", CHARGE_MAP.size());
    }

    /** For testing. */
    public static void main(String[] args)
    {
        load();
    }
} 
