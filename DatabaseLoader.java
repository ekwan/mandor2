import java.util.concurrent.*;
import java.util.concurrent.atomic.*;

/**
 * Provides a quick way to load all the various libraries using multiple threads.
 */
public class DatabaseLoader
{
    private DatabaseLoader()
    {
        throw new IllegalArgumentException("not instantiable");
    }

    public static class LoaderThread extends Thread
    {
        private final AtomicInteger signal;

        public LoaderThread(AtomicInteger signal)
        {
            this.signal = signal;
        }
    }

    /**
     * Loads all libraries with multiple threads.
     */
    public static void go()
    {
        final AtomicInteger signal = new AtomicInteger();
        
        LoaderThread ramachandranThread = new LoaderThread(signal) {
                                        public void run()
                                            {
                                                RamachandranDatabase.load();
                                                signal.getAndIncrement();
                                            } };
        ramachandranThread.start();
            
        LoaderThread rotamerThread = new LoaderThread(signal) {
                                        public void run()
                                            {
                                                RotamerDatabase.load();
                                                signal.getAndIncrement();
                                            } };
        rotamerThread.start();
        
        LoaderThread omegaThread = new LoaderThread(signal) {
                                        public void run()
                                            {
                                                OmegaDatabase.load();
                                                signal.getAndIncrement();
                                            } };
        omegaThread.start();
        
        while (true)
            {
                if ( signal.get() == 3 )
                    break;
                GeneralThreadService.wait(100);
            }
    }
}
