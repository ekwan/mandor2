import java.util.*;
import java.util.concurrent.*;
import java.util.concurrent.atomic.*;
import java.io.*;
import com.google.common.collect.ImmutableList;

/**
 * This singleton runs local work units.
 */
public class GeneralThreadService implements Singleton
{
    /** The thread pool size. */
    public static final int NUMBER_OF_THREADS = Settings.NUMBER_OF_THREADS;

    /** How many jobs can wait in the queue at one time. */
    public static final int JOB_CAPACITY = 10000000;

    private static final CustomThreadPoolExecutor EXECUTOR_SERVICE;

    public static final AtomicInteger numberRunning = new AtomicInteger();

    /** Not instantiable. */
    private GeneralThreadService()
    {
        throw new IllegalArgumentException("not instantiable");
    }

    /** Static initializer. */
    static
    {
        // create a thread pool that has exactly NUMBER_OF_THREADS threads
        // when started, the executor service will create threads up to NUMBER_OF_THREADS
        // these threads will be kept running until the service is shut down
        // if more than JOB_CAPACITY jobs is placed in the queue, then the rejected execution policy will decide what happens
        // the "caller runs policy" returns the excess work to the calling thread

        // threads only created on demand
        EXECUTOR_SERVICE = new CustomThreadPoolExecutor(NUMBER_OF_THREADS, // core pool size
                                                       NUMBER_OF_THREADS, // maximum pool size
                                                       1L, // keep alive time
                                                       TimeUnit.MINUTES, // keep alive time unit
                                                       new ArrayBlockingQueue<Runnable>(JOB_CAPACITY,true), // work queue
                                                       new CustomThreadFactory("thread pool"), // thread factory
                                                       new ThreadPoolExecutor.CallerRunsPolicy()); // rejected execution policy
        System.out.println("GeneralThreadService initialized with " + NUMBER_OF_THREADS + " threads.");
    }

    /**
     * Forces the static initializer to run.
     */
    public static void initialize()
    {
    }

    /**
     * Returns the number of jobs that are waiting in the thread pool queue.
     * @return the number of waiting jobs
     */
    public static int queueSize()
    {
        return EXECUTOR_SERVICE.getQueue().size();
    }

    /**
     * Submits the specified job.
     * @param u the work to run
     * @return a promise for the result in the future
     */
    public static Future<Result> submit(WorkUnit u)
    {
        return EXECUTOR_SERVICE.submit(u);
    }

    /** Alias method. */
    public static Future<RemoteResult> submit(RemoteWorkUnit u)
    {
        return EXECUTOR_SERVICE.submit(u);
    }

    /**
     * The thread pool that will run the work on the local node.
     */
    private static class CustomThreadPoolExecutor extends ThreadPoolExecutor
    {
        public CustomThreadPoolExecutor(int corePoolSize, int maximumPoolSize, long keepAliveTime, TimeUnit unit,
                                        ArrayBlockingQueue<Runnable> workQueue, ThreadFactory threadFactory, RejectedExecutionHandler handler)
        {
            super(corePoolSize, maximumPoolSize, keepAliveTime, unit, workQueue, threadFactory, handler);
        }

        protected void beforeExecute(Thread t, Runnable r)
        {
            super.beforeExecute(t,r);
            numberRunning.getAndIncrement();
        }

        protected void afterExecute(Runnable r, Throwable t)
        {
            super.afterExecute(r,t);
            numberRunning.getAndDecrement();

            // print exceptions if any
            try
                {
                    Future<?> future = (Future<?>)r;
                    future.get();
                }
            catch (ExecutionException e)
                {
                    System.out.println("\n=== EXECUTION EXCEPTION ===\n");
                    e.getCause().printStackTrace();
                }
            catch (Exception e)
                {
                    e.printStackTrace();
                }
        }
    }

    /**
     * Creates the worker threads that will perform the work. 
     */
    private static class CustomThreadFactory implements ThreadFactory
    {
        private final String poolName;

        public CustomThreadFactory(String poolName)
        {
            this.poolName = poolName;
        }

        // instead of creating Threads, create WorkerThreads (a descendent of Thread)
        public Thread newThread(Runnable runnable)
        {
            return new WorkerThread(runnable, poolName);
        }
    }

    /**
     * The worker threads that will actually do the work.
     */
    public static class WorkerThread extends Thread
    {
        private static final AtomicInteger created = new AtomicInteger(); // thread safe integer

        public WorkerThread(Runnable runnable, String name)
        {
            super(runnable, "thread " + created.incrementAndGet());
            setDaemon(true);
        }

        public void run()
        {
            super.run();
        }
    }

    /**
     * Waits for the specified time.
     * @param time time in milliseconds to wait
     */
    public static void wait(int time)
    {
        try
            {
                Thread.sleep(time);
            }
        catch (InterruptedException e)
            {
            }
    }

    /**
     * Waits for the specified jobs to finish.  No progress report.
     * @param futures the futures of the jobs we want to wait for
     */
    public static void silentWaitForFutures(List<Future<Result>> futures)
    {
        int totalJobs = futures.size();
        if ( totalJobs == 0 )
            return;
        while (true)
            {
                int numberDone = 0;
                for (Future<Result> f : futures)
                    if ( f.isDone() )
                        numberDone++;
                if ( numberDone == totalJobs )
                    break;
                wait(50);
            }
    }

    /**
     * Waits for the specified local jobs to finish.  Progress report given.
     */
    public static void waitForFutures(List<Future<Result>> futures)
    {
        int totalJobs = futures.size();
        if ( totalJobs == 0 )
            return;
        while (true)
            {
                int numberDone = 0;
                for (Future<Result> f : futures)
                    if ( f.isDone() )
                        numberDone++;
                System.out.printf("%d of %d jobs complete      queue: %d    \r", numberDone, totalJobs, queueSize());
                if ( numberDone == totalJobs )
                    break;
                wait(50);
            }
    }
}
