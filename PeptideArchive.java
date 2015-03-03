import java.util.*;
import java.io.*;
import com.google.common.collect.*;

public class PeptideArchive implements Serializable
{
    public static final long serialVersionUID = 1L;

    public final Map<String,List<Peptide>> map;

    public PeptideArchive(Map<String,List<Peptide>> map)
    {
        this.map = ImmutableMap.copyOf(map);
    }

    @Override
    public String toString()
    {
        String returnString = "";
        for (String s : map.keySet())
            returnString += String.format("%9s (%d structures)\n", s, map.get(s).size());
        return returnString;
    }

    /**
     * Writes the peptide to disk.  Will not die from exceptions.
     * @param filename the filename to write to
     */
    public void checkpoint(String filename)
    {
        if ( filename == null || filename.length() == 0 )
            throw new NullPointerException("must specify a filename");
        try
            {
                 FileOutputStream fileOut = new FileOutputStream(filename);
                 ObjectOutputStream out = new ObjectOutputStream(fileOut);
                 out.writeObject(this);
                 out.close();
                 fileOut.close();
            }
        catch (Exception e)
            {
                e.printStackTrace();
            }
    }

    /**
     * Loads a peptide checkpoint from disk.  WIll not die from exceptions.
     * @param filename the filename to write to
     * @return the deserialized peptide
     */
    public static PeptideArchive load(String filename)
    {
        if ( filename == null || filename.length() == 0 )
            throw new NullPointerException("must specify a filename");
        try
            {
                FileInputStream fileIn = new FileInputStream(filename);
                ObjectInputStream in = new ObjectInputStream(fileIn);
                PeptideArchive p = (PeptideArchive)in.readObject();
                in.close();
                fileIn.close();
                return p;
            }
        catch (Exception e)
            {
                e.printStackTrace();
                return null;
            }
    }

    public static void main(String[] args)
    {
        PeptideArchive archive = load("checkpoints/designs.chk");
        System.out.println(archive);
    }
}
