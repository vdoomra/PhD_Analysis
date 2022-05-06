#   Author: Vassu Doomra (Stony Brook University) Date: April 17, 2022

#!/apps/bin/python3
from subprocess import call
import sys, os, time, tarfile

def main():

    email = "vdoomra@jlab.org"

    energyscan = [2.2, 4.4, 6.6, 8.8, 11]
    
    for energy in energyscan:

       config = "SlimGeneral_" + str(energy)
       outDir = "/volatile/halla/moller12gev/vdoomra/"+config
       srcDir = "/volatile/halla/moller12gev/vdoomra/moller_optics_"+str(energy)
       if not os.path.exists(outDir):
         os.makedirs(outDir)

       nrStart = 0
       nrStop = 200
       submit = 1

       for jobNr in range(nrStart,nrStop):
           print("Starting job setup for jobID: " + str(jobNr))

           jobName = config + '_%03d'%jobNr
           outDirFull = outDir + "/" + jobName
        
           if not os.path.exists(outDirFull):
             os.makedirs(outDirFull)
    
           createShellScript(srcDir,outDirFull)
           createSBATCHfile(outDirFull,jobName,jobNr)

           if submit==1:
            print("submitting", jobNr)
            call(["sbatch",outDirFull+"/submit.sh"])

    print("All done for config ",config," for #s between ",nrStart, " and ", nrStop)

def createShellScript(srcDir,outDirFull):
    fname = outDirFull + '/slim.sh'
    f=open(fname,"w")
    f.write("#!/bin/sh\n")
    f.write("LIST=`ls -lrt "+ srcDir + "/moller_optics_*/o_remoll.root "+"| "+"awk '{print $9}'`\n")
    f.write("NUM=0\n")
    f.write("for file in $LIST\n")
    f.write("do\n")
    f.write("if (( $NUM == $1 ))\n")
    f.write("then\n")
    f.write("cp /work/halla/moller12gev/vdoomra/optics/remoll/optics_output/SlimGeneral.C .\n")
    f.write("root -l -b << EOF\n")
    f.write(".X SlimGeneral.C(\"$file\",0)\n")
    f.write("EOF\n")
    f.write("fi\n")
    f.write("NUM=$(( $NUM + 1 ))\n")
    f.write("done\n")
    f.close()
    os.system('chmod +x ' + fname)
    return 0

def createSBATCHfile(outDirFull,jobName,jobNr):

    f=open(outDirFull+"/submit.sh","w")
    f.write("#!/bin/bash\n")
    f.write("#SBATCH --ntasks=1\n")
    f.write("#SBATCH --job-name="+jobName+"\n")
    f.write("#SBATCH --output="+outDirFull+"/log.out\n")
    f.write("#SBATCH  --error="+outDirFull+"/log.err\n")
    f.write("#SBATCH --partition=production\n")
    f.write("#SBATCH --account=halla\n")
    f.write("#SBATCH --mem-per-cpu=5000\n")
    f.write("#SBATCH --exclude=farm19104,farm19105,farm19106,farm19107,farm1996,farm19101,farm180131,farm1998,farm1997\n")
    f.write("cd "+outDirFull+"\n")
    f.write("./slim.sh "+str(jobNr)+"\n")
    f.write("rm -rf macros slim.sh SlimGeneral.C\n")
    f.write("rm -r submit.sh\n")
    f.close()
    return 0
  
if __name__ == '__main__':
    main()
