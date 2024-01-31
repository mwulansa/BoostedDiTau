import sys,os
import glob

Sample = sys.argv[1]
filePerList = sys.argv[2]

fileList = glob.glob("sampleList/Ntuple_"+Sample+"*v2024.txt")

outputDir = "filelists/"
filePath = outputDir+Sample

print(fileList)

for f in fileList:
    smallfile = None
    i = 0
    with open(f) as bigfile:
        fileName = "_".join(f.split("/")[1].split("_")[:-1])
        if not os.path.exists(filePath+"/"+fileName):
            os.makedirs(filePath+"/"+fileName)
        for lineno, line, in enumerate(bigfile):
            if lineno % int(filePerList) == 0:
                if smallfile:
                    smallfile.close()
                    i +=1
                small_filename = fileName+"_"+str(int(i))+".txt"                
                smallfile = open(filePath+"/"+fileName+"/"+small_filename, "w")
                print(filePath+"/"+fileName+"/"+small_filename)
            smallfile.write(line)
        if smallfile:
            smallfile.close()
