jobid=$1
#apptsv=$2
#trainOut=$3
#testOut=$4

apptsv="./sampleData/sampled_apps.tsv"
trainOut="./sampleData/sampled_train_apps.tsv"
testOut="./sampleData/sampled_test_apps.tsv"

grep -r "\t$jobid" $apptsv
grep -r "$jobid" ../Dataset/jobs.tsv >> "./sampleData/sampled_jobs.tsv"
grep -r "\t$jobid" $apptsv | grep -r "\tTrain" >> $trainOut
grep -r "\t$jobid" $apptsv | grep -r "\tTest" >> $testOut
