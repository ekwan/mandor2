    cp $1 old_xyz
    awk -f translate.awk $1 > temp.xyz
    mv temp.xyz $1
