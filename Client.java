import java.util.*;
import java.net.*;
import java.io.*;
import java.util.concurrent.*;
import java.util.concurrent.atomic.*;
import com.google.common.collect.*;

/**
 * This class performs RemoteWorkUnits dispatched by a remote Server and sends back RemoteResults.
 */
public class Client implements Singleton
{
    /** Try to connect to the server on this port. */
    public static final int LISTENING_PORT = Server.LISTENING_PORT;

    /** The host name that the Server is running on. */
    public static String HOSTNAME = "enj12.rc.fas.harvard.edu";

    /** Send this string to establish a connection with the Server. */
    public static final String HANDSHAKE = "handshake";

    /** If this string is received, close down the server. */
    public static final String CLOSE = "close";

    ///** How many seconds to wait before declaring the Server unreachable. */
    //public static final int TIMEOUT = 1;

    /** How many times to try connecting to the Server before giving up. */
    public static final int MAX_ATTEMPTS = 5;

    /** How long to wait between connection attempts in seconds. */
    public static final int RETRY_INTERVAL = 2;

    private static Socket connection;
    private static InputStream inputStream;
    private static ObjectInputStream incomingObjectStream;
    private static OutputStream outputStream;
    private static ObjectOutputStream outgoingObjectStream;

    /** Will lock on this before sending anything back to the Server. */
    private static Object sendLock = new Object();

    /** This is a singleton. */
    private Client()
    {
        throw new IllegalArgumentException("you aren't supposed to instantiate this");
    }

    /**
     * Connects to the Server, performs work, and shuts down if requested.  I assume we won't be doing anything on the
     * clients other than running this code.
     */
    public static void main(String[] args)
    {
        GeneralThreadService.initialize();

        // use the first command line argument as the hostname
        if ( args.length > 0 )
            HOSTNAME = args[0];

        int attempts=0;
        while (true)
            {
                try
                    {
                        System.out.print("Attempting to connect to " + HOSTNAME + "...");
                        attempts++;
                        connect();
                        System.out.println("connected to " + connection.getInetAddress() + " and ready to run jobs.");
                        break;
                    }
                catch (InterruptedException e)
                    {
                        System.out.println("Interrupted!");
                    }
                catch (EOFException e)
                    {
                        System.out.println("Broken pipe; reconnecting.");
                    }
                catch (ConnectException | ClassNotFoundException e)
                    {
                        if ( e.getMessage().toLowerCase().indexOf("connection refused") > -1 )
                            System.out.println("connection refused.");
                        else
                            {
                                System.out.println("Error connecting:");
                                e.printStackTrace();
                            }
                    }
                catch (SocketException e)
                    {
                        System.out.println("Unable to connect, retrying.");
                    }
                catch (UnknownHostException e)
                    {
                        System.out.println("unknown host name...aborting!");
                        System.exit(1);
                    }
                catch (Exception e)
                    {
                        e.printStackTrace();
                    }

                // check if we've exceeded the maximum number of connection attempts
                if ( attempts > MAX_ATTEMPTS )
                    {
                        System.out.println("Maximum number of connection attempts exceeded.");
                        System.exit(1);
                        break;
                    }

                // wait before retrying
                System.out.print("Waiting to retry...");
                GeneralThreadService.wait(RETRY_INTERVAL*1000);
                System.out.println("done.\n");

            }

        // start accepting work
        doWork();

        // once doWork exits, we are finished
        System.out.println("Client has shut down.");
    }

    /**
     * Establishes a connection to the server.
     */
    private static void connect() throws IOException, InterruptedException, ConnectException, ClassNotFoundException, SocketException, UnknownHostException, EOFException
    {
        // establish connection
        System.out.print("opening socket...");
        connection = new Socket(HOSTNAME, LISTENING_PORT);
        System.out.println("connection established.");
        //connection.setSoTimeout(TIMEOUT*1000);
        
        // create streams
        outputStream = new BufferedOutputStream(connection.getOutputStream());
        outputStream.flush();
        outgoingObjectStream = new ObjectOutputStream(outputStream);
        outgoingObjectStream.flush();
        inputStream = connection.getInputStream();
        BufferedInputStream tempStream = new BufferedInputStream(inputStream);
        incomingObjectStream = new ObjectInputStream(tempStream);

        // send handshake
        System.out.print("Handshaking...");
        outgoingObjectStream.writeObject(HANDSHAKE);
        outgoingObjectStream.flush();
        System.out.print("handshake sent...");

        // wait for handshake
        Object incomingObject = incomingObjectStream.readObject();
        if ( incomingObject instanceof String )
            {
                String thisHandshake = (String)incomingObject;
                if ( !thisHandshake.equals(HANDSHAKE) )
                    throw new ConnectException("Error handshaking (wrong handshake).");
            }
        else
            throw new ConnectException("Error handshaking (wrong object type).");
    }

    /**
     * Dispatches the work sent by the server.
     */
    private static void doWork()
    {
        while (true)
            {
                try
                    {
                        Object incomingObject = incomingObjectStream.readObject();
                        //System.out.println("object received!");
                        if ( incomingObject instanceof RemoteWorkUnit )
                            {
                                RemoteWorkUnit unit = (RemoteWorkUnit)incomingObject;
                                System.out.printf("Job received (ID %d)\n", unit.getServerID());
                                GeneralThreadService.submit(unit);
                            }
                        else if ( incomingObject instanceof String )
                            {
                                String incomingString = (String)incomingObject;
                                if ( incomingString.equals(CLOSE) )
                                    {
                                        System.out.println("Close request received.");
                                        break;
                                    }
                                else
                                    System.out.println("Unknown string received: " + incomingString);
                            }
                        else if ( incomingObject == null )
                            System.out.println("received a null!");
                        else
                            System.out.println("Unknown object type received: " + incomingObject.getClass());
                    }
                catch (SocketException e)
                    {
                        if ( e.getMessage().indexOf("Connection reset") > -1 )
                            System.out.println("Connection reset.");
                        else
                            e.printStackTrace();
                        break;
                    }
                catch (EOFException e)
                    {
                        break;
                    }
                catch (Exception e)
                    {
                        e.printStackTrace();
                    }
            }

        try
            {
                connection.close();
            }
        catch (IOException e)
            {
                e.printStackTrace();
            }
        System.out.println("Connection to server closed.");
        
        System.exit(0);
    }

    /**
     * Sends the specifed result back to the server.  If we are not on a Client, this has no effect.
     * @param remoteResult what to send back
     */
    public static void sendResult(RemoteResult remoteResult)
    {
        // ensure we are on a client
        if ( !Settings.MAIN_CLASS.equals("Client") )
            {
                System.out.println("Warning, trying to send back a result when not running a Client.  Ignoring.");
                System.out.println(Settings.MAIN_CLASS);
                return;
            }

        // try to send the result back
        int attempts = 0;
        while ( attempts < 5 )
            {
                synchronized(sendLock)
                    {
                        try
                            {
                                attempts++;
                                outgoingObjectStream.writeObject(remoteResult);
                                outgoingObjectStream.flush();
                                outgoingObjectStream.reset();
                                System.out.printf("Sent back result (ID=%d).\n", remoteResult.getServerID());
                                return;
                            }
                        catch (Exception e)
                            {
                                e.printStackTrace();
                                GeneralThreadService.wait(2000);
                            }
                    }
            }
        System.out.println("Giving up on sending back unit.");
    }
}
