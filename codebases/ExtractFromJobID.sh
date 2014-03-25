jobid=$1
apptsv=$2
trainOut=$3
testOut=$4

grep -r "\t$jobid" $apptsv
grep -r "\t$jobid" $apptsv | grep -r "\tTrain" >> $trainOut
grep -r "\t$jobid" $apptsv | grep -r "\tTest" >> $testOut
