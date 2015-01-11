# this script takes one parameter, the three letter abbreviation for the amino acid
#
# this script gives us:
#
# a tinker xyz file containing the amobea atom types, this will go in $1.xyz
# a printout of the oplsaal atom types
#
# it will ensure that the amoeba and oplssal geometries line up
#
# the $1.xyz file can be sent to omnisol to get the surface tensions
# 
# oplsaal is just opls with better parameters for proteins
#
# usage: ./protein.sh ALA

# make AMOEBA template
sed s/RESIDUE/$1/g protein.input > protein.input.temp
sed s/PARAMETERFILE/amoebapro13.prm/g protein.input.temp > protein.input.temp2
protein < protein.input.temp2 > /dev/null
rm protein.input.temp*
rm $1.int
mv $1.xyz $1_amoeba.xyz
rm $1.seq

# make OPLS template
sed s/RESIDUE/$1/g protein.input > protein.input.temp
sed s/PARAMETERFILE/oplsaal.prm/g protein.input.temp > protein.input.temp2
protein < protein.input.temp2 > /dev/null
rm protein.input.temp*
rm $1.int
mv $1.xyz $1_opls.xyz
rm $1.seq

# ensure geometries are identical
awk '{if (NR>1) {print $3, $4, $5}}' $1_amoeba.xyz > geometry_amoeba.txt
awk '{if (NR>1) {print $3, $4, $5}}' $1_opls.xyz > geometry_opls.txt
diff geometry_amoeba.txt geometry_opls.txt > geometry_diff.txt
if [ -s geometry_diff.txt ]; then
    echo geometry mismatch
    exit 1
else
    echo made xyz files, geometries match
fi
rm geometry_amoeba.txt
rm geometry_opls.txt
rm geometry_diff.txt

# get opls atom types
echo "OPLS atom types: atom number, OPLS type"
awk '{if (NR>1) { printf "OPLStype  %4d  %4d\n", $1, $6 }}' $1_opls.xyz

# clean up
mv $1_amoeba.xyz $1.xyz
rm $1_opls.xyz
