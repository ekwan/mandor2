#!/usr/bin/env bash
# usage: benchmark.sh num_of_cores
#
# ./benchmark.sh 6

if [ $# -ne 1 ]; then
    echo "Must specify number of cores."
    exit 1
fi

hostname=`hostname | awk 'BEGIN{FS="."}{print $1}'`
echo Hostname is $hostname.
rm -f $hostname*

date
start=`date +%s`

# submit jobs

for i in $(seq 1 $1); do
    cp original.xyz ${hostname}_$i.xyz
    #./minimize ${hostname}_$i.xyz amoebapro13.prm 0.2 > ${hostname}_$i.out &
    #echo ./analyze ${hostname}_$i.xyz amoebapro13.prm d ${hostname}_$i.out 
    analyze ${hostname}_$i.xyz amoebapro13.prm d > ${hostname}_$i.out &
done

# wait for jobs to complete

wait

#while :; do
#    counter=0
#    for i in ${hostname}_*.out; do
#        #grep -q "Normal Termination" $i &> /dev/null
#        grep -q "Total Potential Energy" $i &> /dev/null
#        if [ "$?" == "0" ]; then
#            let "counter=counter+1"
#        fi
#    done
#    echo -ne $counter of $1" jobs complete\r"
#    if [ "$counter" == "$1" ]; then
#        break;
#    fi
#
#    sleep 0.01
#done

echo
end=`date +%s`
date
runtime=$((end-start))
echo $runtime seconds
