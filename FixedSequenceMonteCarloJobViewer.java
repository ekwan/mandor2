import java.util.*;
import java.io.*;

public class FixedSequenceMonteCarloJobViewer
{
    public static FixedSequenceMonteCarloJob load(String filename)
    {
        if ( filename == null || filename.length() == 0 )
            throw new NullPointerException("must specify a filename");
        try
            {
                FileInputStream fileIn = new FileInputStream(filename);
                ObjectInputStream in = new ObjectInputStream(fileIn);
                Object object = in.readObject();
                if ( ! ( object instanceof FixedSequenceMonteCarloJob ) )
                    return null;
                FixedSequenceMonteCarloJob job = (FixedSequenceMonteCarloJob)object;
                in.close();
                fileIn.close();
                return job;
            }
        catch (Exception e)
            {
                e.printStackTrace();
            } 
        return null;
    }

    public static void main(String args[])
    {
        int index = 0;
        for (File f : new File("checkpoints").listFiles())
            {
                String filename = f.getName();
                if ( filename.endsWith(".chk") )
                    {
                        System.out.printf("Loading %s...", filename);
                        FixedSequenceMonteCarloJob job = load("checkpoints/" + filename);
                        if ( job == null )
                            {
                                System.out.println("error.");
                                continue;
                            }
                        List<Peptide> list = job.bestPeptides.getList();
                        String prefix = String.format("checkpoints/dump_%02d_", index);
                        Peptide.writeGJFs(list, prefix, 2, 5);
                        System.out.println("done.");
                        index++;
                    }
            }
    }
}
