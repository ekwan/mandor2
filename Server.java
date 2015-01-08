import java.io.*;
import java.net.*;
import java.util.*;
import java.util.concurrent.*;
import com.google.common.collect.*;
import java.util.concurrent.atomic.*;

/**
 * This is a general purpose class that send work units to Clients.  If connections are lost,
 * the units will be automatically requeued.  Submit a job with Server.submit(RemoteWorkUnit).
 * This will return a RemoteFuture<RemoteResult>.  This result can be down-cast to the appropriate
 * type by the caller.
 */
public class Server implements Singleton
{
    // Networking Parameters

    /** Will listen to connecitons on this port. */
    public static final int LISTENING_PORT = 32007;

    /** Seconds to wait before declaring a computer unreachable. */
    public static final int TIMEOUT = 1;

    /** How many jobs to assign per client.  Maps hostname substrings to number of jobs. */
    public static final Map<String,Integer> HOSTNAME_DATABASE;

    /** How many jobs to assign per client by default. */
    public static final int DEFAULT_JOBS_PER_HOST = 2;

    /** If set to true, will override the specified capacity for the ConnectionThreads. */
    public static boolean SHOULD_OVERRIDE_CAPACITY = false;

    /** Each ConnectionThread will run this number of jobs instead. */
    public static final int OVERRIDE_CAPACITY = 8;

    /** All the live connections and how many jobs have been sent out to each one. */
    private static final Map<ConnectionThread,Integer> CONNECTIONS;

    /** Keeps track of the state of all the running work units. */
    private static final WorkUnitDatabase DATABASE;

    /** Keeps track of the work units that have not been run yet. */
    private static Map<RemoteWorkUnit,RemoteFuture<RemoteResult>> QUEUE;

    /** The thread that keeps track of the incoming connections and parcels out the connections. */
    private static final MasterThread MASTER_THREAD;

    /** Static initializer. */
    static
    {
        CONNECTIONS = Collections.synchronizedMap(new HashMap<ConnectionThread,Integer>());
        DATABASE = new WorkUnitDatabase();
        QUEUE = Collections.synchronizedMap(new LinkedHashMap<RemoteWorkUnit,RemoteFuture<RemoteResult>>());

        // Add to this map to change the number of jobs per client.  Hostnames matching
        // the substrings stored in the keys will be set to have a maximum number of jobs
        // equalling the integer in the value.  If no value is set, the default value of
        // DEFAULT_JOBS_PER_HOST will be used.
        Map<String,Integer> tempMap = new HashMap<>();
        tempMap.put("dae", 2);
        tempMap.put("enj", 3);    
        HOSTNAME_DATABASE = ImmutableMap.copyOf(tempMap);
        
        MASTER_THREAD = new MasterThread();
        MASTER_THREAD.setDaemon(true);
    }

    private Server()
    {
        throw new IllegalArgumentException("not instantiable");
    }

    /**
     * Start listening for connections.
     */
    public static void initialize()
    {
        MASTER_THREAD.start();
    }

    /**
     * Dispatches work to the ConnectionThreads and requeues work.
     */
    private static class MasterThread extends Thread
    {
        public MasterThread()
        {
        }

        public void run()
        {
            // start accepting connections and wait for all jobs to complete
            ServerSocket listener = null;
            Socket connection = null;
            ConnectionThread connectionThread = null;
            
            System.out.println("Listening on port " + LISTENING_PORT + "...");

            while ( true )
                {
                    // check for dead threads and requeue any dead units
                    List<ConnectionThread> deadThreads = new ArrayList<>();
                    synchronized(CONNECTIONS)
                        {
                            for (ConnectionThread thread : CONNECTIONS.keySet())
                                {
                                    if ( thread.isAlive() )
                                        continue;
                                    deadThreads.add(thread);
                                    Map<RemoteWorkUnit,RemoteFuture<RemoteResult>> deadUnits = DATABASE.markAsDead(thread);
                                    if ( deadUnits.size() > 0 )
                                        {
                                            Map<RemoteWorkUnit,RemoteFuture<RemoteResult>> tempQueue = QUEUE;
                                            QUEUE = Collections.synchronizedMap(new LinkedHashMap<RemoteWorkUnit,RemoteFuture<RemoteResult>>());
                                            QUEUE.putAll(deadUnits);
                                            QUEUE.putAll(tempQueue);
                                            System.out.printf("Connection to %s lost, so requeued %d units.\n", thread.address, deadUnits.size());
                                        }
                                }
                        }
                    for (ConnectionThread t : deadThreads)
                        CONNECTIONS.remove(t);

                    // dispatch new work
                    int count = 0;
                    while ( QUEUE.size() > 0 )
                        {
                            count++;
                            // see if there is someone free to do some work
                            ConnectionThread thread = getNextFreeThread();
                            if ( thread == null )
                                {
                                    //System.out.println("There is work to do, but there are no free threads to do it with.");
                                    break;
                                }

                            // get the next piece of work
                            RemoteWorkUnit unit = null;
                            RemoteFuture<RemoteResult> future = null;
                            synchronized (QUEUE)
                            {
                                for (RemoteWorkUnit u : QUEUE.keySet())
                                    {
                                        unit = u;
                                        future = QUEUE.get(u);
                                        break;
                                    }
                            }
                            if ( unit == null || future == null )
                                throw new NullPointerException("unexpected null");

                            // dispatch it and mark it as started
                            try
                                {
                                    thread.dispatch(unit);
                                }
                            catch (Exception e)
                                {
                                    e.printStackTrace();
                                    continue;
                                }
                            DATABASE.markAsStarted(unit,thread,future);
                            QUEUE.remove(unit);
                            int jobsRunning = CONNECTIONS.get(thread);
                            jobsRunning++;
                            CONNECTIONS.put(thread,jobsRunning);

                            // don't get stuck here
                            if (count > 10)
                                break;
                        }

                    // listen for new connections
                    try
                        {
                            listener = new ServerSocket(LISTENING_PORT);
                            listener.setSoTimeout(TIMEOUT*1000);
                            connection = listener.accept();
                            String hostname = connection.getInetAddress().getCanonicalHostName();
                            
                            System.out.printf("[ %s ] Opened a socket to %s (%s).\n", new Date().toString(), hostname, connection.getInetAddress());
                            listener.close();
                            
                            // set the maximum number of jobs each host can handle
                            int jobCapacity = DEFAULT_JOBS_PER_HOST;
                            for (String s : HOSTNAME_DATABASE.keySet())
                                {
                                    if ( hostname.indexOf(s) > -1 )
                                        {
                                            jobCapacity = HOSTNAME_DATABASE.get(s);
                                            break;
                                        }
                                }
                            connectionThread = new ConnectionThread(connection, jobCapacity);
                            connectionThread.start();
                        }
                    catch (BindException e)
                        {
                            if (e.getMessage().equals("Address already in use"))
                                System.out.println("A server is already running on this port!");
                            else
                                e.printStackTrace();
                        }
                    catch (SocketTimeoutException e)
                        {
                            try
                                {
                                    if ( listener != null )
                                        listener.close();
                                }
                            catch (Exception e2)
                                {
                                    e2.printStackTrace();
                                    System.exit(1);
                                }
                        }
                    catch (ConnectException e)
                        {
                            System.out.println(e.getMessage());
                        }
                    catch (EOFException e)
                        {
                            System.out.println("Connection to " + connection.getInetAddress() + " closed unexpectedly.");
                            CONNECTIONS.remove(connectionThread);
                        }
                    catch (Exception e)
                        {
                            e.printStackTrace();
                            break;
                        }
                }
            System.out.println("Server has shut down.");
        }
    }

    /**
     * Represents a connection between the Server and a Client.
     */
    private static class ConnectionThread extends Thread
    {
        private final Socket connection;
        private InputStream incomingStream;
        private ObjectInputStream incomingObjectStream;
        private OutputStream outgoingStream;
        private ObjectOutputStream outgoingObjectStream;

        /** The host name. */
        private final String address;

        /** If set to true, will terminate at the next opportunity. */
        private final AtomicBoolean closeSignal;

        /** How many jobs this thread can run at once. */
        private final int jobCapacity;

        /** To avoid concurrency issues with sending objects. */
        private final Object sendLock = new Object();

        /** For establishing connections. */
        public static final String HANDSHAKE = "handshake";

        /** For closing connections. */
        public static final String CLOSE = "close";
    
        /**
         * Tries to make a connection with the specified socket.
         * @param connection the socket to connect to
         */
        public ConnectionThread(Socket connection, int jobCapacity) throws IOException, ConnectException, EOFException, ClassNotFoundException
        {
            // set streams
            this.connection = connection;
            outgoingStream = new BufferedOutputStream(connection.getOutputStream());
            outgoingStream.flush();
            outgoingObjectStream = new ObjectOutputStream(outgoingStream);
            outgoingObjectStream.flush();
            incomingStream = new BufferedInputStream(connection.getInputStream());
            incomingObjectStream = new ObjectInputStream(incomingStream);

            // receive handshake
            Object incomingObject = incomingObjectStream.readObject();
            if ( incomingObject instanceof String )
                {
                    String thisString = (String)incomingObject;
                    if ( ! thisString.equals(HANDSHAKE) )
                        throw new ConnectException("Error handshaking -- wrong text!");
                }
            else
                throw new ConnectException("Error handshaking -- wrong object type!");
        
            // send handshake
            outgoingObjectStream.writeObject(HANDSHAKE);
            outgoingObjectStream.flush();

            // get the host name
            address = connection.getInetAddress().getCanonicalHostName();
            
            // mark this as a live connection, set the number of running jobs on this thread to zero
            CONNECTIONS.put(this,0);
            this.jobCapacity = jobCapacity;
            
            // when set to true, will close connection
            closeSignal = new AtomicBoolean(false);
            System.out.printf("Connection to %s established with a maximum of %d jobs.\n", address, jobCapacity);
        }

        /**
         * Waits for incoming objects.  If they are results, the pile of results are updated. 
         */
        public void run()
        {
            while (true)
                {
                    try
                        {
                            Object incomingObject = incomingObjectStream.readObject();
                            if ( incomingObject instanceof RemoteResult )
                                {
                                    // cast to RemoteResult
                                    RemoteResult remoteResult = (RemoteResult)incomingObject;
                                    
                                    // mark as finished
                                    Long ID = Long.valueOf(remoteResult.getServerID());
                                    System.out.printf("Received work unit ID %d.  %d units remain queued.\n", ID, QUEUE.size());
                                    WorkUnitDatabase.DatabaseEntry databaseEntry = DATABASE.markAsFinished(ID);
                                    ConnectionThread thread = databaseEntry.thread;
                                    synchronized (CONNECTIONS)
                                        {
                                            int jobsRunning = CONNECTIONS.get(thread);
                                            jobsRunning--;
                                            CONNECTIONS.put(thread,jobsRunning);
                                        }

                                    // set the result of this job 
                                    RemoteFuture<RemoteResult> futureTask = databaseEntry.remoteFuture;
                                    futureTask.set(remoteResult);
                                }
                            else
                                System.out.printf("Unknown object type received from %s.\n", address);
                        }
                    catch (SocketTimeoutException e)
                        {
                            if ( closeSignal.get() )
                                {
                                    System.out.printf("Connection to %s closed gracefully.\n", address);
                                    break;
                                }
                        }
                    catch (EOFException | SocketException e)
                        {
                            System.out.printf("Connection to %s closed due to exception %s.\n", address, e.getMessage());
                            break;
                        }
                    catch (Exception e)
                        {
                            e.printStackTrace();
                            break;
                        }
                }
        }

        /**
         * Sends the specified work unit to the remote host.
         */
        public void dispatch(RemoteWorkUnit unit) throws IOException
        {
            synchronized(sendLock)
                {
                    outgoingObjectStream.writeObject(unit);
                    outgoingObjectStream.flush();
                    outgoingObjectStream.reset();
                }
        }

        /**
         * Returns the internet address of this connection.
         * @return the host name
         */
        public String getHostName()
        {
            return address;
        }

        /**
         * Gracefully closes this conneciton at the next opportunity.
         */
        public void closeConnection()
        {
            closeSignal.set(true);
        }
    }

    /**
     * Keeps track of which work units are currently running on other machines.  When a work unit is started,
     * call @link{#markAsStarted(RemoteWorkUnit,ConnectionThread,RemoteFuture)}.  When a work unit has finished, call
     * this method: @link{#markAsFinished(Long)} .  Work units that die while they are running due to a connection failure
     * will be requeued by Server.
     */
    private static class WorkUnitDatabase
    {
        /** Maps IDs to work units and the threads that are running them. */
        private final Map<Long,DatabaseEntry> map;

        /** Synchronizes the parallel lists. */
        private final Object internalLock;

        /** Standard constructor. */
        public WorkUnitDatabase()
        {
            map = new HashMap<Long,DatabaseEntry>();
            internalLock = new Object();
        }

        /** Represents a piece of work that is in progress on a remote computer. */
        private static class DatabaseEntry
        {
            public final RemoteWorkUnit unit;
            public final ConnectionThread thread;
            public final RemoteFuture<RemoteResult> remoteFuture;

            public DatabaseEntry(RemoteWorkUnit unit, ConnectionThread thread, RemoteFuture<RemoteResult> remoteFuture)
            {
                this.unit = unit;
                this.thread = thread;
                this.remoteFuture = remoteFuture;
            }

            @Override
            public String toString()
            {
                return String.format("%s (%s)", unit.toString(), thread.address);
            }

            @Override
            public int hashCode()
            {
                return Objects.hash(unit,thread);
            }

            @Override
            public boolean equals(Object obj)
            {
                if ( obj == null )
                    return false;
                if ( obj == this ) 
                    return true;
                if ( !(obj instanceof DatabaseEntry) )
                    return false;

                DatabaseEntry d = (DatabaseEntry)obj;
                if ( unit.equals(d.unit) && thread.equals(d.thread) )
                    return true;
                return false;
            }
        }

        /**
         * Mark the specified work unit as being in progress.
         * @param unit the work unit
         * @param thread the thread that is running the work unit
         * @param futureTask where to put the result when the calculation is complete
         * @return information about what is running the job
         */
        public DatabaseEntry markAsStarted(RemoteWorkUnit unit, ConnectionThread thread, RemoteFuture<RemoteResult> futureTask)
        {
            if ( unit == null || thread == null )
                throw new NullPointerException("markAsStarted does not allow nulls!");
            synchronized (internalLock)
                {
                    Long ID = Long.valueOf(unit.getServerID());
                    if ( map.containsKey(ID) )
                        throw new IllegalArgumentException("duplicate server ID key");
                    DatabaseEntry newEntry = new DatabaseEntry(unit,thread,futureTask);
                    map.put(ID, newEntry);
                    return newEntry;
                }
        }

        /**
         * Mark the specified work unit as finished by removing it from the database.
         * @param ID the ID of the completed work unit
         * @return information about what was running the job
         */
        public DatabaseEntry markAsFinished(Long ID)
        {
            DatabaseEntry entry = null;
            synchronized (internalLock)
                {
                    entry = map.remove(ID);
                }
            if ( entry == null )
                throw new NullPointerException("no match found for ID " + ID);
            return entry;
        }

        /**
         * Marks all the work by a specified thread as dead.  Returns the dead units and their futures so they can be re-queued.
         * @param thread the thread we are interested in
         * @return a map from the units to their futures
         */
        public Map<RemoteWorkUnit,RemoteFuture<RemoteResult>> markAsDead(ConnectionThread thread)
        {
            Map<RemoteWorkUnit,RemoteFuture<RemoteResult>> returnMap = new HashMap<>();
            List<DatabaseEntry> deadEntries = new ArrayList<>();
            synchronized(internalLock)
                {
                    for (DatabaseEntry entry : map.values())
                        {
                            RemoteWorkUnit unit = entry.unit;
                            ConnectionThread t = entry.thread;
                            RemoteFuture<RemoteResult> remoteFuture = entry.remoteFuture;
                            if ( t == thread )
                                {
                                    returnMap.put(unit,remoteFuture);
                                    deadEntries.add(entry);
                                }
                        }
                    for (DatabaseEntry entry : deadEntries)
                        map.values().remove(entry);
                }
            return returnMap;
        }

        /**
         * Returns the number of jobs running remotely.
         * @return the number of jobs
         */
        public int size()
        {
            return map.size();
        }
    }

    /**
     * Runs the requested work unit.  The unit is actually placed in a queue.  When a thread becomes free to do work,
     * the unit will be sent to a client.  When the result comes back, the result of the future will be set.
     * @param unit the unit to run
     * @return a ticket promising the result at a later time
     */
    public static RemoteFuture<RemoteResult> submit(RemoteWorkUnit unit)
    {
        RemoteFuture<RemoteResult> future = new RemoteFuture<>();
        QUEUE.put(unit, future);
        return future;
    }

    /**
     * Finds a free ConnectionThread to do some work.  If none are available, returns null.
     * @return the free ConnectionThread, if it exists.
     */
    public static ConnectionThread getNextFreeThread()
    {
        synchronized(CONNECTIONS)
            {
                for (ConnectionThread t : CONNECTIONS.keySet())
                    {
                        Integer jobCapacity = t.jobCapacity;
                        if ( SHOULD_OVERRIDE_CAPACITY )
                            jobCapacity = OVERRIDE_CAPACITY;
                        Integer jobsRunning = CONNECTIONS.get(t);
                        if ( jobCapacity > jobsRunning )
                            return t;
                    }
            }
        return null;
    }

    /**
     * Waits for the specified remote jobs to finish.  Progress report given. 
     */
    public static void waitForFutures(List<RemoteFuture<RemoteResult>> futures)
    {
        int totalJobs = futures.size();
        String lastString = "";
        while (true)
            {
                int numberDone = 0;
                for (Future<RemoteResult> f : futures)
                    if ( f.isDone() )
                        numberDone++;
                int dispatched = DATABASE.size(); 
                int queueSize = QUEUE.size();
                String reportString = String.format("%d of %d jobs complete      dispatched: %d    queued: %d    \r", numberDone, totalJobs, dispatched, queueSize);
                if ( !reportString.equals(lastString) )
                    {
                        System.out.printf(reportString);
                        lastString = reportString;
                    }
                if ( numberDone == totalJobs )
                    break;
                GeneralThreadService.wait(250);
            }
    }
}    
