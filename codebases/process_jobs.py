import csv
import re
with open('./sampleData/sampled_jobs.tsv','rb') as tsvin, open('./sampleData/newJobs.csv', 'wb') as csvout:
    tsvin = csv.reader(tsvin, delimiter='\t')
    csvout = csv.writer(csvout)

    # process header
    header = tsvin.next()
    nColumns = len(header)

    # preparation for non-textual discretization
    dictionaries = [] # dictionary for discretizing each column
    values = [] # numerical count of distinct label
    for i in range(0, nColumns):
        dictionaries.append({})
        values.append(0)

    # preparation for textual discretization
    keywords = [r'manage', r'market', r'bilingual',
                r'computer', r'design', r'technician', r'Business',
                r'bachelor', r'Administrative', r'programer|programming']
    nKeywords = len(keywords)

    row = header 
    csvout.writerows([[row[0]]+ [row[2]] + row[5:9]] + row[11:])

    discretSet = [2] + range(5,9) # set of index
    desIndex = 3
    for row in tsvin:
        # non-textual feature
        for i in discretSet:
            if not dictionaries[i].has_key(row[i]):
                dictionaries[i][row[i]] = values[i]
                values[i] += 1
            row[i] = dictionaries[i][row[i]]
        # textual feature
        description = row[desIndex]
        tfeat = []
        for i in range(0, nKeywords):
            result = re.search(keywords[i], description)
            if result is not None:
                tfeat.append(1)
            else:
                tfeat.append(0)
        # write to out file
        jobFeature = [row[0]]+ [row[2]] + row[5:9] + tfeat
        print jobFeature
        csvout.writerows([jobFeature])

