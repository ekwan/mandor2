# For technical reasons, amino acids are not allowed to have any duplicate atoms.
# Duplicate atoms are defined as any with the same symbol and position.
# This program shifts the coordinates in a file by a random amount.
# usage: awk -f translate.awk input.xyz > output.xyz
BEGIN {
    randomShift = 2.0
    srand()
    }
{
    if (NR == 1)
        {
            randomX = rand()*randomShift
            randomY = rand()*randomShift
            randomZ = rand()*randomShift
            print
        }
    else if (length($0) > 0)
        {
            x = $3 + randomX
            y = $4 + randomY
            z = $5 + randomZ
            printf("%6d  %-6s%10.6f  %10.6f  %10.6f%6d", $1, $2, x, y, z, $6);
            if (NF>6)
                {
                    for (i=7; i <= NF; i++)
                        printf("%6d", $i)
                }
            printf("\n")
        }
    else
        {
            print
        }
}
