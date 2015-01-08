# This awk script splits the Dunbrack Ramachadran database data into separate files.
# Usage: awk -f split_NDRD.awk NDRD_whatever.txt
# Output will be dumped to "NDRD_split....txt".  You should delete any pre-existing
{
    if ( $1 == "#" )
        outputFilename = "NDRD_split_header.txt"
    else
        outputFilename = "NDRD_split_" $1 ".txt"
    print $0 >> outputFilename
    if ( NR % 100 == 0 )
        printf "%20s %d\r", outputFilename, NR
}
END {printf "\n"}
