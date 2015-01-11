# converts a tinker xyz file to omnisol input file geometry block
{
    if ( NR > 1 )
        {
            symbol = substr($2,1,1)
            printf "%-2s %10.6f %10.6f %10.6f\n", symbol, $3, $4, $5
        }
}
