import java.util.*;
import java.io.*;
import com.google.common.collect.*;
import org.apache.commons.math3.geometry.euclidean.threed.*;
import org.jgrapht.*;
import org.jgrapht.graph.*;
import org.jgrapht.alg.*;
import java.util.concurrent.*;

/**
 * This class creates peptides.
 */
public class PeptideFactory 
{
    /** This class is not instantiable. */
    private PeptideFactory()
    {
        throw new IllegalArgumentException("not instantiable");
    }

    /**
     * Sets hairpin angles for the given peptide.
     * This method expects a single D-proline followed by a glycine.
     * @param peptide the peptide to adjust the angles for
     * @return the adjusted peptide
     */
    public static Peptide setHairpinAngles(Peptide peptide)
    {
        // figure out which residues to adjust
        Residue prolineResidue = null;
        Residue glycineResidue = null;
        Integer prolineResidueIndex = null;
        Integer glycineResidueIndex = null;
        for (int i=0; i < peptide.sequence.size(); i++)
            {
                Residue residue = peptide.sequence.get(i);
                if ( residue.aminoAcid == AminoAcid.DPRO )
                    {
                        if ( prolineResidue != null )
                            throw new IllegalArgumentException("found two D-prolines");
                        prolineResidue = residue;
                        glycineResidue = peptide.sequence.get(i+1);
                        prolineResidueIndex = i;
                        glycineResidueIndex = i+1;
                        if ( glycineResidue.aminoAcid != AminoAcid.GLY )
                            throw new IllegalArgumentException("glycine expected after D-proline");
                        i++;
                    }
            }
        if ( prolineResidue == null || glycineResidue == null )
            throw new NullPointerException("hairpin not found");

        // set the dihedral angles
        // these angles are from the Balaram crystal structure
        Peptide newPeptide = BackboneMutator.setOmega (peptide, prolineResidueIndex,  175.0);
        newPeptide         = BackboneMutator.setPhiPsi(peptide, prolineResidueIndex,   53.0, -132.0);
        newPeptide         = RotamerMutator.setChis   (peptide, prolineResidueIndex,  ImmutableList.of(-22.0, 32.0, -29.0));
        newPeptide         = BackboneMutator.setOmega (peptide, glycineResidueIndex, -178.0);
        newPeptide         = BackboneMutator.setPhiPsi(peptide, glycineResidueIndex,  -96.0,    9.0);
        
        // set hairpin flags
        List<Residue> newSequence = new ArrayList<>(newPeptide.sequence.size());
        for (int i=0; i < newPeptide.sequence.size(); i++)
            {
                Residue residue = newPeptide.sequence.get(i);
                if ( i == prolineResidueIndex || i == glycineResidueIndex )
                    newSequence.add(residue.setHairpin(true));
                else
                    newSequence.add(residue.setHairpin(false));
            }
        newPeptide = new Peptide(newPeptide.name, newPeptide.contents, newPeptide.connectivity, newSequence, newPeptide.energyBreakdown);
        return newPeptide;
    }

    /**
     * Builds a peptide given a sequence of ProtoAminoAcids.
     * @param inputSequence the requested sequence specified in the N to C direction
     * @return the Peptide that embodies the geometry and metadata of the sequence
     */
    public static Peptide createPeptide(List<ProtoAminoAcid> inputSequence)
    {
        if ( inputSequence.size() < 2 )
            throw new IllegalArgumentException("A peptide must have at least two residues in it.");
        
        // create temporary lists for creating a new Residue
        List<AminoAcid>                                 aminoAcids              = new LinkedList<>();
        List<ProtoTorsion>                              omegas                  = new LinkedList<>();
        List<ProtoTorsion>                              phis                    = new LinkedList<>();
        List<ProtoTorsion>                              psis                    = new LinkedList<>();
        List<List<ProtoTorsion>>                        chis                    = new LinkedList<>();
        List<Atom>                                      HNs                     = new LinkedList<>();
        List<Atom>                                      Ns                      = new LinkedList<>();
        List<Atom>                                      Os                      = new LinkedList<>();
        List<Atom>                                      Cs                      = new LinkedList<>();
        List<Atom>                                      CAs                     = new LinkedList<>();
        List<Atom>                                      HAs                     = new LinkedList<>();
        List<String>                                    descriptions            = new LinkedList<>();
        List<Pair<Atom,Atom>>                           prochiralConnections    = new LinkedList<>(); 
        List<List<Atom>>                                atoms                   = new LinkedList<>();
        List<Boolean>                                   isHairpins              = new LinkedList<>();

        // for Molecule part of new Peptide
        String                                          newName                 = "";
        List<Atom>                                      newContents             = new LinkedList<>();
        SimpleWeightedGraph<Atom,DefaultWeightedEdge>   newConnectivity         = new SimpleWeightedGraph<>(DefaultWeightedEdge.class);
        Set<Pair<Atom,Atom>>                            newBonds                = new HashSet<>();

        // keep track of connections so new bonds can be formed
        List<Pair<Atom,Atom>>                           NStickyConnections      = new LinkedList<>();
        List<Pair<Atom,Atom>>                           CStickyConnections      = new LinkedList<>();
        List<Pair<Atom,Atom>>                           newConnections          = new LinkedList<>();
        List<Pair<Integer,Integer>>                     newIndexConnections     = new LinkedList<>();

        // build peptide in N to C direction
        for (int i=0; i < inputSequence.size(); i++)
            {
                // get data for this amino acid
                ProtoAminoAcid p = inputSequence.get(i).shift(i+1);

                AminoAcid                               tempAminoAcid           = p.residue.aminoAcid;
                ProtoTorsion                            tempOmega               = p.residue.omega;
                ProtoTorsion                            tempPhi                 = p.residue.phi;
                ProtoTorsion                            tempPsi                 = p.residue.psi;
                List<ProtoTorsion>                      tempChis                = p.residue.chis;
                Atom                                    tempHN                  = p.residue.HN;
                Atom                                    tempN                   = p.residue.N;
                Atom                                    tempO                   = p.residue.O;
                Atom                                    tempC                   = p.residue.C;
                Atom                                    tempCA                  = p.residue.CA;
                Atom                                    tempHA                  = p.residue.HA;
                String                                  tempDescription         = p.residue.description;
                Pair<Atom,Atom>                         tempNStickyConnection   = p.NStickyConnection;
                Pair<Atom,Atom>                         tempCStickyConnection   = p.CStickyConnection;
                Pair<Atom,Atom>                         tempProchiralConnection = p.residue.prochiralConnection;
                List<Atom>                              tempAtoms               = p.residue.atoms;
                boolean                                 tempIsHairpin           = false;

                Molecule                                tempMolecule            = p.molecule;

                // append to name
                if ( i == 0 )
                    newName += String.format("%5s", tempAminoAcid.shortName);
                else
                    newName += String.format(" - %5s", tempAminoAcid.shortName);
                
                aminoAcids.add(tempAminoAcid);

                // add relevant atoms from this ProtoAminoAcid
                tempAtoms = new LinkedList<>(tempAtoms);
                if ( i > 0 )
                    {
                        // delete N-terminal cap
                        Atom keepAtom = tempNStickyConnection.getFirst();
                        Atom deleteAtom = tempNStickyConnection.getSecond();
                        Set<Atom> toBeDeleted = tempMolecule.getHalfGraph(keepAtom, deleteAtom);
                        tempAtoms.removeAll(toBeDeleted);
                    }

                if ( i < inputSequence.size() - 1 )
                    {
                        // delete C-terminal cap
                        Atom keepAtom = tempCStickyConnection.getFirst();
                        Atom deleteAtom = tempCStickyConnection.getSecond();
                        Set<Atom> toBeDeleted = tempMolecule.getHalfGraph(keepAtom, deleteAtom);
                        tempAtoms.removeAll(toBeDeleted);
                    }
                newContents.addAll(tempAtoms);  // the atoms in the whole peptide
                atoms.add(tempAtoms);           // list of the atoms in each residue

                // add relevant connectivity from each ProtoAminoAcid
                for (DefaultWeightedEdge e : tempMolecule.connectivity.edgeSet())
                    {
                        // assumes all edge weights are 1
                        Atom fromAtom = tempMolecule.connectivity.getEdgeSource(e);
                        Atom toAtom   = tempMolecule.connectivity.getEdgeTarget(e);
                        if ( tempAtoms.contains(fromAtom) && tempAtoms.contains(toAtom) )
                            newBonds.add( new Pair<Atom,Atom>(fromAtom, toAtom) );
                    }

                // form new bond
                NStickyConnections.add(tempNStickyConnection);
                CStickyConnections.add(tempCStickyConnection);

                if ( i > 0 )  // no need to form a bond on the first residue
                    {
                        // join the C sticky atom of the previous residue
                        // to the N sticky atom of the current residue
                        // sticky atom connections are defined as ordered pairs: (keep, delete)
                        Atom lastCStickyAtom = CStickyConnections.get(i-1).getFirst();
                        Atom thisNStickyAtom = NStickyConnections.get(i).getFirst();
                        Pair<Atom,Atom> newConnection = new Pair<>(lastCStickyAtom, thisNStickyAtom);
                        //System.out.println(newContents.indexOf(lastCStickyAtom)+1);
                        //System.out.println(newContents.indexOf(thisNStickyAtom)+1);

                        // mark this in a list of all the new peptide bonds
                        newConnections.add(newConnection);

                        // mark for addition to connectivity graph
                        newBonds.add(newConnection);

                        // specify the new connections by atom number because atoms will be moved
                        // indices: 0, 1, 2, ..., N
                        Integer fromAtomIndex = newContents.indexOf(lastCStickyAtom);
                        Integer toAtomIndex = newContents.indexOf(thisNStickyAtom);
                        Pair<Integer,Integer> newIndexPair = new Pair<>(fromAtomIndex, toAtomIndex);
                        newIndexConnections.add(newIndexPair);
                    }

                // update backbone angles
                if ( i == 0 )
                    {
                        // copy ProtoTorsions for first residue
                        omegas.add(tempOmega);
                        phis.add(tempPhi);
                        psis.add(tempPsi);
                    }
                else if ( i > 0 )
                    {
                        // get last ProtoTorsions
                        ProtoTorsion lastPsi = psis.remove(psis.size()-1);
                        
                        // update ProtoTorsions for previous residue
                        ProtoTorsion newLastPsi = new ProtoTorsion(lastPsi.atom1, lastPsi.atom2, lastPsi.atom3, tempOmega.atom3);
                        psis.add(newLastPsi);

                        // update ProtoTorsions for current residue
                        ProtoTorsion newOmega = new ProtoTorsion(lastPsi.atom2, lastPsi.atom3, tempOmega.atom3, tempOmega.atom4);
                        ProtoTorsion newPhi   = new ProtoTorsion(lastPsi.atom3, tempPhi.atom2,   tempPhi.atom3,   tempPhi.atom4);
                        omegas.add(newOmega);
                        phis.add(newPhi);
                        psis.add(tempPsi);
                    }

                // add sidechain torsion angles
                chis.add(tempChis);

                // update special atom fields
                HNs.add(tempHN); 
                Ns.add(tempN);
                Os.add(tempO);
                Cs.add(tempC);
                CAs.add(tempCA);
                HAs.add(tempHA);

                // add other fields
                descriptions.add(tempDescription);
                prochiralConnections.add(tempProchiralConnection);
                isHairpins.add(tempIsHairpin);
            }

        // create new connectivity graph
        for (Atom a : newContents)
            newConnectivity.addVertex(a);
        for (Pair<Atom,Atom> p : newBonds)
            {
                // assumes all edges have a weight of 1.0
                Atom fromAtom = p.getFirst();
                Atom toAtom = p.getSecond();
                DefaultWeightedEdge e = newConnectivity.addEdge(fromAtom, toAtom);
                newConnectivity.setEdgeWeight(e, 1.0);
                //System.out.println( (newContents.indexOf(fromAtom)+1) + " -- " + (newContents.indexOf(toAtom)+1) );
            }

        // create residue objects to make a Sequence
        List<Residue> sequence = new LinkedList<>();
        for (int i=0; i < inputSequence.size(); i++)
            {
                Residue residue = new Residue(aminoAcids.get(i), omegas.get(i), phis.get(i), psis.get(i), chis.get(i),
                                              HNs.get(i), Ns.get(i), Os.get(i), Cs.get(i), CAs.get(i), HAs.get(i),
                                              descriptions.get(i), prochiralConnections.get(i), atoms.get(i), isHairpins.get(i));

                sequence.add(residue);
            }
        sequence = ImmutableList.copyOf(sequence);

        // create new peptide object
        Peptide peptide = new Peptide(newName, newContents, newConnectivity, sequence, EnergyBreakdown.BLANK);

        // adjust distances
        for (Pair<Integer,Integer> bond : newIndexConnections)
            {
                Atom fromAtom = peptide.contents.get(bond.getFirst());
                Atom toAtom   = peptide.contents.get(bond.getSecond());
                peptide  = peptide.setMolecule( peptide.setDistance(fromAtom,toAtom,1.32) );
            }
        
        // set peptide bonds to sp2
        for (int i=0; i < inputSequence.size(); i++)
            {
                // set amide carbonyl carbon to sp2
                Residue r = peptide.sequence.get(i);
                ProtoTorsion omega = r.omega;
                Atom fromAtom = omega.atom2;
                Atom toAtom = omega.atom3;
                peptide = peptide.setMolecule( peptide.set_sp2(fromAtom,toAtom) );

                // set amide nitrogen to sp2
                r = peptide.sequence.get(i);
                omega = r.omega;
                fromAtom = omega.atom3;
                toAtom = omega.atom2;

                // only force (amide substituent 1 -- N -- amide substituent 2) angle to 120
                // if this is not a proline
                if ( r.aminoAcid.isProline() )
                    peptide = peptide.setMolecule( peptide.set_sp2(fromAtom,toAtom,false) );
                else
                    peptide = peptide.setMolecule( peptide.set_sp2(fromAtom,toAtom, true) );
            }

        // set all backbone torsions to trans
        for (int i=0; i < inputSequence.size(); i++)
            peptide = BackboneMutator.setOmega(peptide, i, 180.0);

        return peptide;
    }

    /** For testing. */
    public static void main(String[] args)
    {
        DatabaseLoader.go();
        List<ProtoAminoAcid> sequence = ProtoAminoAcidDatabase.getSpecificSequence("arg","met","standard_ala","gly","d_proline", "gly", "phe", "val", "hd", "l_pro");
        Peptide peptide = PeptideFactory.createPeptide(sequence);
        
        for (int i=0; i < peptide.sequence.size(); i++)
            {
                peptide = BackboneMutator.mutateOmega(peptide, i);
                peptide = BackboneMutator.mutatePhiPsi(peptide, i);
                peptide = RotamerMutator.mutateChis(peptide, i);
            }
        peptide = PeptideFactory.setHairpinAngles(peptide);

        GaussianInputFile f = new GaussianInputFile(peptide);
        f.write("test_peptides/test.gjf");

        //peptide = BackboneMutator.mutateOmega(peptide, 1);
        //peptide = BackboneMutator.mutatePhiPsi(peptide, peptide.sequence.get(1));
        //f = new GaussianInputFile(peptide);
        //f.write("test_peptides/test2.gjf");
    }
}  
