import sys, os

#Sample = sys.argv[1]
#filePerList = sys.argv[2]
#process = sys.argv[3]

Sample = sys.argv[1]
filePerList = sys.argv[2]
#process = sys.argv[3]

smallfile = None
i = 0

def checkAndMakeDir(dir):
    if not os.path.exists(dir):
        os.mkdir(dir)

process = Sample.split("_")

if process[0] == "DYJetsToLL":
    sampleName = process[0]+"_"+process[1]
    part = process[2]
    outputDir = "filelists/"+sampleName+"/"+part+"/"
    sfname = sampleName+"_"+part
elif process[0] == "QCD" or process[0] == "WJetsToLNu":
    sampleName = process[0]
    part = process[1]
    outputDir = "filelists/"+sampleName+"/"+part+"/"
    sfname = sampleName+"_"+part
elif process[0] == "ST":
    if process[2] == "4f": 
        sampleName = process[0]+"_"+process[1]
    else: sampleName = process[0]+"_"+process[1]+"_"+process[2]
    outputDir = "filelists/"+sampleName+"/"
    sfname = sampleName
else:
    sampleName = process[0]
    outputDir = "filelists/"+sampleName+"/"
    sfname = sampleName


# sampleName = "TTToSemiLeptonic"
# sfname = sampleName

# outputDir = "filelists/"+sampleName+"/"

print(outputDir)
checkAndMakeDir("filelists/"+sampleName)
checkAndMakeDir(outputDir)

#with open("sampleList/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8.txt") as bigfile:

with open("sampleList/"+Sample+".txt") as bigfile:
    for lineno, line, in enumerate(bigfile):
        if lineno % int(filePerList) == 0:
            if smallfile:
                smallfile.close()
                i +=1
            small_filename = sfname+"_"+str(int(i))+".txt"
            smallfile = open(outputDir+small_filename, "w")
            print(outputDir+small_filename)
        smallfile.write(line)
    if smallfile:
        smallfile.close()
