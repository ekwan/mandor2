# creates a beta hairpin with the angles in hairpin.input and puts it in hairpin.pdb
rm -f hairpin2.int
rm -f hairpin2.xyz
rm -f hairpin2.seq
rm -f hairpin2.pdb
protein < hairpin2.input
xyzpdb hairpin2.xyz amoebapro13.prm
rm -f hairpin2.int
rm -f hairpin2.seq
