# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

baseFileName = "../neura-sparse01/edit-v3.csv"
suppFileName = "../neura-sparse01/edit-v3dynamic.csv"
outFileName = "../neura-sparse01/edit-v3merge.csv"

baseF = open(baseFileName, "r")
suppF = open(suppFileName, "r")
outF = open(outFileName, "w")

suppF.readline()
line2 = suppF.readline()
words2 = line2.split(",")
words2[1] = words2[1].replace("-imuStepDetect.csv", "").replace("-revStepDetect.csv", "")

for line in baseF:
    words = line.split(",")
    words[1] = words[1].replace("-imuStepDetect.csv", "").replace("-revStepDetect.csv", "")
    if words[0] == "save" and words2[1] == words[1]:
        # append data from words[1] supplement
        print(words2[1])
        while True:
            line2 = suppF.readline()
            if line2 == "": break
            words2 = line2.split(",")
            words2[1] = words2[1].replace("-imuStepDetect.csv", "").replace("-revStepDetect.csv", "")
            if words2[0] == "load": 
                continue
            elif words2[0] == "save":
                line2 = suppF.readline()
                if line2 != "":
                    words2 = line2.split(",")
                    words2[1] = words2[1].replace("-imuStepDetect.csv", "").replace("-revStepDetect.csv", "")
                break
            else:
                outF.write(line2)
    outF.write(line)
    
baseF.close()
suppF.close()
outF.close()