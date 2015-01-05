# creates a beta hairpin with the angles in hairpin.input and puts it in hairpin.pdb
rm -f hairpin.int
rm -f hairpin.xyz
rm -f hairpin.seq
rm -f hairpin.pdb
protein < hairpin.input
xyzpdb hairpin.xyz amoebapro13.prm
rm -f hairpin.int
rm -f hairpin.xyz
rm -f hairpin.seq
