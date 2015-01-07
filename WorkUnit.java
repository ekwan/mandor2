import java.util.concurrent.*;

/**
 * This call represents a piece of work that can be done in parallel.
 */
public interface WorkUnit extends Callable<Result>
{
    public Result call();
}
