import csv

with open('./sampleData/sampled_train_apps.tsv') as trainIn, \
        open('./sampleData/sampled_test_apps.tsv') as testIn, \
        open('./sampleData/newUsers.csv') as userIn, \
        open('./sampleData/newJobs.csv') as jobIn, \
        open('./sampleData/joint_train.csv', 'w') as csvtrainout, \
        open('./sampleData/joint_test.csv', 'w') as csvtestout:

    csvtrainout = csv.writer(csvtrainout, delimiter=' ')
    csvtestout = csv.writer(csvtestout, delimiter=' ')

    trainSet = csv.reader(trainIn, delimiter='\t')
    testSet = csv.reader(testIn, delimiter='\t')
    userSet = csv.reader(userIn, delimiter=',')
    jobSet = csv.reader(jobIn, delimiter=',')

    nTrainColumns = len(trainSet.next())
    nTestColumns = len(testSet.next())
    nUserColumns = len(userSet.next())
    nJobColumns = len(jobSet.next())

    # set up user dictionary
    allusers = {}
    for user in userSet:
        ukey = int(user[0])
        uvalue = tuple(user[3:])
        # print ukey, uvalue 
        allusers.update({ukey:uvalue})
        
    # set up job dictionary
    alljobs = {}
    for job in jobSet:
        jkey = int(job[0])
        jvalue = tuple(job[1:])
        # print jkey, jvalue
        alljobs.update({jkey:jvalue})

    # set up training join table
    trainEntities = []
    trainLabels = []
    for app in trainSet:
        userid = int(app[0])
        if allusers.has_key(userid):
            userfeat = list(allusers[userid])
        else:
            print "ERROR"

        jobid = int(app[-1])
        if alljobs.has_key(jobid):
            jobfeat = list(alljobs[jobid])
        else:
            print "ERROR"

        jointEntity = ["+1"]
        #jointEntity = [jointEntity] + [userid, jobid]
        for i in range(len(userfeat)):
            if userfeat[i] is None or len(userfeat[i]) == 0:
                userfeat[i] = 0
            jointEntity.append(str(i+1) + ":" + str(userfeat[i]))
        nUserFeatures = len(userfeat)
        for j in range(len(jobfeat)):
            if jobfeat[j] is None or len(jobfeat[j]) == 0:
                jobfeat[j] = 0
            jointEntity.append(str(nUserFeatures+j+1) + ":" + str(jobfeat[j]))
        #print jointEntity
        #trainEntities.append(jointEntity)
        #trainLabels.append(1)
        csvtrainout.writerows([jointEntity])

    # set up test join table
    testEntities = []
    testLabels = []
    for app in testSet:
        userid = int(app[0])
        if allusers.has_key(userid):
            userfeat = list(allusers[userid])
        else:
            print "ERROR"

        jobid = int(app[-1])
        if alljobs.has_key(jobid):
            jobfeat = list(alljobs[jobid])
        else:
            print "ERROR"

        jointEntity = ["+1"]
        #jointEntity = [jointEntity] + [userid, jobid]
        for i in range(len(userfeat)):
            if userfeat[i] is None or len(userfeat[i]) == 0:
                userfeat[i] = 0
            jointEntity.append(str(i+1) + ":" + str(userfeat[i]))
        nUserFeatures = len(userfeat)
        for j in range(len(jobfeat)):
            if jobfeat[j] is None or len(jobfeat[j]) == 0:
                jobfeat[j] = 0
            jointEntity.append(str(nUserFeatures+j+1) + ":" + str(jobfeat[j]))
        #print jointEntity
        #testEntities.append(jointEntity)
        #testLabels.append(1)
        csvtestout.writerows([jointEntity])

