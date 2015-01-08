Neighbor-dependent Ramachandran Distributions 
Version 1: released April 29, 2010

D. Ting, M. Shapovalov, G. Guoli, R. Mitra, Michael I. Jordan, and R.L.Dunbrack, Jr.
Fox Chase Cancer Center
Philadelphia PA 19111, USA
215 728 2434
http://dunbrack.fccc.edu/ndrd
Roland.Dunbrack@fccc.edu

User must have license and accepted its terms when
obtaining and using these data. The files can not be redistributed to
other individuals outside your own lab group without permission. If
you want to include these maps in software for further distribution,
then please contact us at Roland.Dunbrack@fccc.edu. It can be
worked out in a reasonably easy manner.

Each file contains probabilities of phi,psi for specific residue
types given the residue type of a neighbor to the left or to the
right. "ALL" means all neighbor residues of the residue in question
were kept in the calculation. The format is this:


So for instance, here are some lines for LEU-right-PRO Ramachandran distribution

Res Dir   Neigh  phi  psi  Probability   log(prob) Cumulative_sum
LEU right PRO  -175  -130  5.237406e-08   16.76485  2.082167e-04
LEU right PRO  -175  -125  4.726758e-08   16.86744  2.082640e-04
LEU right PRO  -175  -120  4.449988e-08   16.92778  2.083085e-04
LEU right PRO  -175  -115  4.332896e-08   16.95444  2.083518e-04
LEU right PRO  -175  -110  4.208621e-08   16.98355  2.083939e-04
LEU right PRO  -175  -105  3.898562e-08   17.06007  2.084329e-04
LEU right PRO  -175  -100  3.459067e-08   17.17968  2.084675e-04
LEU right PRO  -175   -95  3.050099e-08   17.30551  2.084980e-04
LEU right PRO  -175   -90  2.744314e-08   17.41115  2.085254e-04
LEU right PRO  -175   -85  2.587084e-08   17.47015  2.085513e-04
LEU right PRO  -175   -80  2.586474e-08   17.47039  2.085771e-04

Res = the residue type for the Ramachandran Distribution
Dir = the direction of the neighbor ("ALL" is ALL residue types at once)
phi and psi = the centers of 5x5 regions
Probability = the probability in the 5x5 regions centered that phi,psi
log(prob) = followed by the log probaility.
Cumulative_sum  = the cumulative sum and can be used for drawing random values from the probability distributions. The sum is 1.0 for each neighbor map.
CPR is cis proline as a central amino acid

To calculate probabilies for triplets, use:

log p(C,L,R) = log p(C,L) + log p(C,R) + const
Once log p(C,L,R) is calculated, calculate p(C,L,R) = exp(log(p(C,L,R)))
Then sum them up for each Ramchandan map, and normalize the probabilities by dividing by the sum.

There are four distribution files:

NDRD_TCBIG.txt = data from Turn, Coil, Bridge, PiHelix, and 310 Helix
NDRD_TCB.txt = data from Turn, Coil and Bridge
NDRD_Conly.txt = data from Coil only (note: small data set)
NDRD_Tonly.txt = data from Turn only (note: small data set)

