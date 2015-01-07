import java.util.*;
import java.util.concurrent.*;

/**
 * This singleton contains all of the backbone-dependent rotamer library data we have.
 */
public class RotamerDatabase implements Singleton
{
    /** This map stores the rotamer libraries for the normal amino acids, which always have a trans omega dihedral angle. */
    private static final Map<AminoAcid,RotamerLibrary> DATABASE;

    /** The rotamer library for cis proline. */
    private static final RotamerLibrary CIS_PROLINE_LIBRARY;

    /** The rotamer lirary for trans proline. */
    private static final RotamerLibrary TRANS_PROLINE_LIBRARY;

    // static initializer
    static
    {
        // make temporary map
        Map<AminoAcid, SidechainRotamerLibrary> tempMap = new ConcurrentHashMap<>();
        System.out.println("Loading rotamer libraries...");

        // load libraries in parallel
        List<Future<Result>> futures = new ArrayList<>();
        AtomicReference<RotamerLibrary> cisTemp = new AtomicReference<>();
        AtomicReference<RotamerLibrary> transTemp = new AtomicReference();
        for (AminoAcid a : AminoAcid.values())
            {
                String filename = null;
                if ( a == AminoAcid.LPRO )
                    {
                        String filename = Settings.ROTAMER_LIBRARY_DIRECTORY + "tpr.bbdep.rotamers.lib";
                        Future<Result> f = loadSpecialLibrary(transTemp, filename);
                        futures.add(f);
                        filename = Settings.ROTAMER_LIBRARY_DIRECTORY + "cpr.bbdep.rotamers.lib";
                        f = loadSpecialLibrary(cisTemp, filename);
                        futures.add(f);
                    }
                else if ( a.rotamerType == AminoAcid.RotamerType.HAS_NO_ROTAMERS ||
                          a.rotamerType == AminoAcid.RotamerType.SPECIAL            )
                    {
                        // do nothing
                    }
                else if ( a.rotamerType == AminoAcid.RotamerType.IS_ROTAMERIC )
                    {
                        String filename = Settings.ROTAMER_LIBRARY_DIRECTORY + a.shortName.toLowerCase() + ".bbdep.rotamers.lib"
                        Future<Result> f = loadLibrary(tempMap, filename);
                        futures.add(f);
                    }
                else if ( a.rotamerType == AminoAcid.RotamerType.NON_ROTAMERIC )
                    {
                        String filename = Settings.ROTAMER_LIBRARY_DIRECTORY + a.shortName.toLowerCase() + ".bbdep.densities.lib";
                        Future<Result> f = loadLibrary(tempMap, filename);
                        futures.add(f);
                    }
                else
                    throw new IllegalArgumentException("unreachable");
            }
        GeneralThreadService.silentWaitForFutures(futures);

        // return immutable copy
        if ( cisTemp.get() == null )
            throw new NullPointerException("unable to load cis proline library");
        if ( transTemp.get() == null )
            throw new NullPointerException("unable to load trans proline library");
        CIS_PROLINE_LIBRARY = cisTemp.get();
        TRANS_PROLINE_LIBRARY = transTemp.get();
        MAP = ImmutableMap.copyOf(tempMap);

        System.out.printf("%d databases have been loaded.\n", DATABASE.size() + 2);
    }

    /** Not instantiable. */
    private RotamerDatabase()
    {
        throw new IllegalArgumentException("not instantiable");
    }

    private static Future<Result> loadLibrary(Map<AminoAcid,RotamerLibrary> map, String filename)
    {
    }

    private static Future<Result> loadSpecialLibrary(AtomicReference<RotamerLibrary> reference, String filename)
    {
    }

    public static List<List<Double>> getRotamerSpace(AminoAcid aminoAcid, double omega, double phi, double psi)
    {
    }

    public static List<List<Double>> getRotamerSpace(AminoAcid aminoAcid, double phi, double psi)
    {
    }

    public static List<Double> getRandomRotamer(AminoAcid aminoAcid, double omega, double phi, double psi)
    {
    }

    public static List<Double> getRandomRotamer(AminoAcid aminoAcid, double phi, double psi)
    {
    }

}
