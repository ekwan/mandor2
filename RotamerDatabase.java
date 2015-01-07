import java.util.*;
import java.util.concurrent.*;
import java.util.concurrent.atomic.*;
import com.google.common.collect.*;

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
        // load libraries in parallel
        System.out.println("Loading rotamer libraries...");
        Map<AminoAcid, RotamerLibrary> tempMap = new ConcurrentHashMap<>();
        List<Future<Result>> futures = new ArrayList<>();
        AtomicReference<RotamerLibrary> cisTemp = new AtomicReference<>();
        AtomicReference<RotamerLibrary> transTemp = new AtomicReference<>();
        for (AminoAcid a : AminoAcid.values())
            {
                if ( a == AminoAcid.LPRO )
                    {
                        String filename = Settings.ROTAMER_LIBRARY_DIRECTORY + "tpr.bbdep.rotamers.lib";
                        SpecialUnit unit = new SpecialUnit(a, filename, transTemp);
                        Future<Result> f = GeneralThreadService.submit(unit);
                        futures.add(f);
                        filename = Settings.ROTAMER_LIBRARY_DIRECTORY + "cpr.bbdep.rotamers.lib";
                        unit = new SpecialUnit(a, filename, cisTemp);
                        f = GeneralThreadService.submit(unit);
                        futures.add(f);
                    }
                else if ( a.rotamerType == AminoAcid.RotamerType.HAS_NO_ROTAMERS ||
                          a.rotamerType == AminoAcid.RotamerType.SPECIAL            )
                    {
                        // do nothing
                    }
                else if ( a.rotamerType == AminoAcid.RotamerType.IS_ROTAMERIC )
                    {
                        String filename = Settings.ROTAMER_LIBRARY_DIRECTORY + a.shortName.toLowerCase() + ".bbdep.rotamers.lib";
                        RegularUnit unit = new RegularUnit(a, filename, tempMap); 
                        Future<Result> f = GeneralThreadService.submit(unit);
                        futures.add(f);
                    }
                else if ( a.rotamerType == AminoAcid.RotamerType.NON_ROTAMERIC )
                    {
                        String filename = Settings.ROTAMER_LIBRARY_DIRECTORY + a.shortName.toLowerCase() + ".bbdep.densities.lib";
                        RegularUnit unit = new RegularUnit(a, filename, tempMap);
                        Future<Result> f = GeneralThreadService.submit(unit);
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
        DATABASE = ImmutableMap.copyOf(tempMap);

        System.out.printf("%d databases have been loaded.\n", DATABASE.size() + 2);
    }

    /** Not instantiable. */
    private RotamerDatabase()
    {
        throw new IllegalArgumentException("not instantiable");
    }

    public static List<List<Double>> getRotamerSpace(AminoAcid aminoAcid, double omega, double phi, double psi)
    {
        return null;
    }

    public static List<List<Double>> getRotamerSpace(AminoAcid aminoAcid, double phi, double psi)
    {
        return null;
    }

    public static List<Double> getRandomRotamer(AminoAcid aminoAcid, double omega, double phi, double psi)
    {
        return null;
    }

    public static List<Double> getRandomRotamer(AminoAcid aminoAcid, double phi, double psi)
    {
        return null;
    }

    private static class RegularUnit implements WorkUnit
    {
        public final AminoAcid aminoAcid;
        public final String filename;
        public final Map<AminoAcid,RotamerLibrary> targetMap;

        public RegularUnit(AminoAcid aminoAcid, String filename, Map<AminoAcid,RotamerLibrary> targetMap)
        {
            this.aminoAcid = aminoAcid;
            this.filename = filename;
            this.targetMap = targetMap;
        }

        public Result call()
        {
            if ( targetMap.containsKey(aminoAcid) )
                throw new IllegalArgumentException("already loaded library for " + aminoAcid.shortName);
            RotamerLibrary library = null;
            if ( aminoAcid.rotamerType == AminoAcid.RotamerType.IS_ROTAMERIC )
                library = new RotamericLibrary(aminoAcid, filename);
            else if ( aminoAcid.rotamerType == AminoAcid.RotamerType.NON_ROTAMERIC )
                library = new NonRotamericLibrary(aminoAcid, filename);
            targetMap.put(aminoAcid, library);

            if ( library == null )
                throw new NullPointerException("library is null");
            targetMap.put(aminoAcid, library);
            return null;
        }
    }

    private static class SpecialUnit implements WorkUnit
    {
        public final AminoAcid aminoAcid;
        public final String filename;
        public final AtomicReference<RotamerLibrary> targetReference;
        
        public SpecialUnit(AminoAcid aminoAcid, String filename, AtomicReference<RotamerLibrary> targetReference)
        {
            this.aminoAcid = aminoAcid;
            this.filename = filename;
            this.targetReference = targetReference;
        }

        public Result call()
        {
            if ( targetReference.get() != null )
                throw new IllegalArgumentException("reference already set");
            RotamerLibrary library = new RotamericLibrary(aminoAcid, filename);
            targetReference.set(library);
            return null;
        }
    }

    public static void load()
    {
        System.out.printf("Database loaded with %d entries.\n", DATABASE.size()+2);
        System.out.println(CIS_PROLINE_LIBRARY);
        System.out.println(TRANS_PROLINE_LIBRARY);
        System.out.println(DATABASE.keySet());
    }

    public static void main(String[] args)
    {
        RotamerDatabase.load();
    }
}
