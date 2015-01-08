import java.util.concurrent.*;

/**
 * Utility class for implementing Future<V>.
 */
public class RemoteFuture<V> implements Future<V>
{
    private V result;

    public RemoteFuture()
    {
        result = null;
    }

    public boolean cancel(boolean mayInterruptIfRunning)
    {
        throw new IllegalArgumentException("operation not supported");
    }

    public V get()
    {
        if ( result == null )
            throw new NullPointerException("result not set");
        return result;
    }

    public V get(long timeout, TimeUnit unit)
    {
        throw new IllegalArgumentException("operation not supported");
    }

    public boolean resultIsNull()
    {
        return result == null;
    }

    public void forget()
    {
        result = null;
    }

    public boolean isCancelled()
    {
        return false;
    }

    public boolean isDone()
    {
        if ( result == null )
            return false;
        return true;
    }

    public void set(V result)
    {
        if ( this.result == null )
            this.result = result;
        else
            throw new IllegalArgumentException("result already set");
    }
}
