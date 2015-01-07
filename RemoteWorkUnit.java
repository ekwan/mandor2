import java.io.*;
import java.util.concurrent.*;
import java.util.concurrent.atomic.*;

/**
 * A work unit designed to be run remotely.
 */
public interface RemoteWorkUnit extends Callable<RemoteResult>, Serializable
{
    /** For serialization. */
    public static final long serialVersionUID = 1L;

    /** For assigning unique IDs. */
    public static final AtomicLong ID_GENERATOR = new AtomicLong();

    /** Performs the work. This should send the result back by using {@link Client#sendResult(RemoteResult)}. */
    public RemoteResult call();

    /** Provides a text description of this work unit. */
    public String toString();

    /**
     * Provides a unique server-generated number that identifies this work unit so that its
     * corresponding Result can be put into the right pile.
     * @return the serial number of this work unit
     */
    public long getServerID();
}
