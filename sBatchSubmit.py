#-----------------------------------------------#
 # Written for the Optics Study                  #
  # Author: Vassu Doomra, Stony Brook University  #
   #                                               #
    #-----------------------------------------------#

#!/apps/bin/python3
from subprocess import call
import sys, os, time, tarfile

def main():

    email = "vdoomra@jlab.org"

    sourceDir = "/work/halla/moller12gev/vdoomra/optics/remoll"
   
    activeDetectors = [28, 32, 33, 34, 35, 270, 9000, 9001, 9002, 9003, 9004, 9005, 
                       9006, 9007, 9008, 9009, 9010, 9011, 9012, 9013, 9014,
                       9015]
    energyscan = [2.2, 4.4, 6.6, 8.8, 11]

    make_tarfile(sourceDir)
   
    for energy in energyscan:

       config = "moller_optics_" + str(energy)
       print("creating " + config)
       outDir = "/volatile/halla/moller12gev/vdoomra/"+config
       if not os.path.exists(outDir):
        os.makedirs(outDir)

       nrEv   = 10000 #100000
       nrStart= 0
       nrStop = 200 #(nrStop -nrStart) 3000
       submit  = 1 ## submit is 1 to submit job, submit is 0 to create folder without submitting the jobs

       print('Running ' + str(nrEv*(nrStop - nrStart)) + ' events...')

       jobName=config + '_%03dkEv'%(nrEv/1000)
       call(["cp",sourceDir+"/jobs/default.tar.gz",
          outDir+"/default.tar.gz"])

       for jobNr in range(nrStart,nrStop): # repeat for jobNr jobs
           print("Starting job setup for jobID: " + str(jobNr))

           jobFullName = jobName + '_%04d'%jobNr
           outDirFull=outDir+"/"+jobFullName
           createMacFile(sourceDir, outDirFull,nrEv,jobNr,activeDetectors,energy)

           createSBATCHfile(sourceDir,outDirFull,jobName,jobNr)

           if submit==1:
            print("submitting", jobName)
            call(["sbatch",sourceDir+"/jobs/"+jobName+".sh"])

    print("All done for config ",config," for #s between ",nrStart, " and ", nrStop)


def createMacFile(srcDir, outDirFull,nrEv,jobNr, detectorList,energy):

    if not os.path.exists(outDirFull):
        os.makedirs(outDirFull)

    f=open(outDirFull+"/"+"macro.mac",'w')
    f.write("/remoll/setgeofile geometry/mollerMother.gdml\n")
    f.write("/remoll/physlist/register QGSP_BERT_HP\n")
    f.write("/remoll/physlist/parallel/enable\n")
    f.write("/remoll/parallel/setfile geometry/mollerParallel.gdml\n")
    f.write("/run/initialize\n")
    f.write("/remoll/addfield "+srcDir+"/map_directory/V2U.1a.50cm.parallel.txt\n")
    f.write("/remoll/addfield "+srcDir+"/map_directory/V2DSg.9.75cm.parallel.txt\n")
    f.write("/remoll/evgen/set elasticC12\n")
    f.write("/remoll/rasx 0 mm\n")
    f.write("/remoll/rasy 0 mm\n")
    f.write("/remoll/evgen/thmin 0.1 deg\n")
    f.write("/remoll/evgen/thmax 2.0 deg\n")
  
    f.write("/remoll/beamene "+ str(energy)+" GeV\n")
    f.write("/remoll/beamcurr 1 microampere\n")
    f.write("/control/execute macros/target/Optics1.mac\n")
    f.write("/control/execute macros/sieve/sieve_in.mac\n")
    f.write("/remoll/SD/disable_all\n")

    for det in detectorList:
        f.write("/remoll/SD/enable "+str(det)+"\n")
        f.write("/remoll/SD/detect lowenergyneutral "+str(det)+"\n")
        f.write("/remoll/SD/detect secondaries "+str(det)+"\n")

    f.write("/remoll/kryptonite/enable\n")
    f.write("/process/list\n")
    f.write("/remoll/seed "+str(int(time.clock_gettime(0)) + jobNr)+"\n")
    f.write("/remoll/filename o_remoll.root\n")
    f.write("/run/beamOn "+str(nrEv)+"\n")
    f.close()
    return 0

def createSBATCHfile(sourceDir,outDirFull,jobName,jobNr):

    if not os.path.exists(sourceDir+"/jobs"):
        os.makedirs(sourceDir+"/jobs")

    f=open(sourceDir+"/jobs/"+jobName+".sh","w")
    f.write("#!/bin/bash\n")
    f.write("#SBATCH --ntasks=1\n")
    f.write("#SBATCH --job-name="+jobName+'_%03d'%jobNr+"\n")
    f.write("#SBATCH --output="+outDirFull+"/log.out\n")
    f.write("#SBATCH  --error="+outDirFull+"/log.err\n")
    f.write("#SBATCH --partition=production\n")
    f.write("#SBATCH --account=halla\n")
    f.write("#SBATCH --mem-per-cpu=5000\n")
    f.write("#SBATCH --exclude=farm19104,farm19105,farm19106,farm19107,farm1996,farm19101\n")
    f.write("cd "+outDirFull+"\n")
    f.write("cp ../default.tar.gz ./\n")
    f.write("tar -zxvf default.tar.gz\n")
    f.write("./remoll macro.mac\n")
    f.write("rm -rf default.tar.gz geometry libremoll.so macro.mac remoll\n")
    f.close()
    return 0

def make_tarfile(sourceDir):
    print("making geometry tarball")
    if os.path.isfile(sourceDir+"/jobs/default.tar.gz"):
        os.remove(sourceDir+"/jobs/default.tar.gz")
    tar = tarfile.open(sourceDir+"/jobs/default.tar.gz","w:gz")
    tar.add(sourceDir+"/build/remoll",arcname="remoll")
    tar.add(sourceDir+"/build/libremoll.so",arcname="libremoll.so")
    tar.add(sourceDir+"/macros/target/Optics1.mac", arcname="macros/target/Optics1.mac")
    tar.add(sourceDir+"/macros/target/Optics2.mac", arcname="macros/target/Optics2.mac")
    tar.add(sourceDir+"/macros/sieve/sieve_in.mac", arcname="macros/sieve/sieve_in.mac")
    tar.add(sourceDir+"/geometry/materials.xml",arcname="geometry/materials.xml")
    tar.add(sourceDir+"/geometry/matrices.xml",arcname="geometry/matrices.xml")
    tar.add(sourceDir+"/geometry/positions.xml",arcname="geometry/positions.xml")
    tar.add(sourceDir+"/geometry/solids/world.xml",arcname="geometry/solids/world.xml")
    tar.add(sourceDir+"/geometry/mollerParallel.gdml" ,arcname="geometry/mollerParallel.gdml") 
    tar.add(sourceDir+"/geometry/mollerMother.gdml" ,arcname="geometry/mollerMother.gdml") 
    tar.add(sourceDir+"/geometry/target/subTargetRegion.gdml" ,arcname="geometry/target/subTargetRegion.gdml") 
    tar.add(sourceDir+"/geometry/target/targetLadder.gdml" ,arcname="geometry/target/targetLadder.gdml") 
    tar.add(sourceDir+"/geometry/electronics/subSBSbunker.gdml" ,arcname="geometry/electronics/subSBSbunker.gdml")
    tar.add(sourceDir+"/geometry/electronics/subPSbunker.gdml" ,arcname="geometry/electronics/subPSbunker.gdml")
    tar.add(sourceDir+"/geometry/hall/hallDaughter_dump.gdml" ,arcname="geometry/hall/hallDaughter_dump.gdml")
    tar.add(sourceDir+"/geometry/hall/subDumpDiffuser.gdml" ,arcname="geometry/hall/subDumpDiffuser.gdml")
    tar.add(sourceDir+"/geometry/upstream/sieve.gdml" ,arcname="geometry/upstream/sieve.gdml")
    tar.add(sourceDir+"/geometry/upstream/upstreamDaughter_merged.gdml" ,arcname="geometry/upstream/upstreamDaughter_merged.gdml")
    tar.add(sourceDir+"/geometry/upstream/upstream_nose_shield_beampipe.gdml" ,arcname="geometry/upstream/upstream_nose_shield_beampipe.gdml")
    tar.add(sourceDir+"/geometry/upstream/inner_upstream_nose_shield_beampipe.gdml" ,arcname="geometry/upstream/inner_upstream_nose_shield_beampipe.gdml")
    tar.add(sourceDir+"/geometry/upstream/upstream.gdml" ,arcname="geometry/upstream/upstream.gdml")
    tar.add(sourceDir+"/geometry/upstream/upstreamToroid.gdml" ,arcname="geometry/upstream/upstreamToroid.gdml")
    tar.add(sourceDir+"/geometry/upstream/upstreamTorusRegion.gdml" ,arcname="geometry/upstream/upstreamTorusRegion.gdml")
    tar.add(sourceDir+"/geometry/upstream/upstream_Conf8.gdml" ,arcname="geometry/upstream/upstream_Conf8.gdml")
    tar.add(sourceDir+"/geometry/upstream/upstreamBeampipe.gdml" ,arcname="geometry/upstream/upstreamBeampipe.gdml")
    tar.add(sourceDir+"/geometry/hybrid/hybridToroid.gdml" ,arcname="geometry/hybrid/hybridToroid.gdml")
    tar.add(sourceDir+"/geometry/hybrid/hybridToroidSupport.gdml" ,arcname="geometry/hybrid/hybridToroidSupport.gdml")
    tar.add(sourceDir+"/geometry/hybrid/hybridDaughter_merged.gdml" ,arcname="geometry/hybrid/hybridDaughter_merged.gdml")
    tar.add(sourceDir+"/geometry/hybrid/hybridDaughter_unified.gdml" ,arcname="geometry/hybrid/hybridDaughter_unified.gdml")
    tar.add(sourceDir+"/geometry/huts/lefthut.gdml" ,arcname="geometry/huts/lefthut.gdml")
    tar.add(sourceDir+"/geometry/showermax/showerMaxGen.gdml" ,arcname="geometry/showermax/showerMaxGen.gdml")
    tar.add(sourceDir+"/geometry/donut/donutConcreteLead.gdml" ,arcname="geometry/donut/donutConcreteLead.gdml")
    tar.add(sourceDir+"/geometry/pion/pionDetectorSystem.gdml" ,arcname="geometry/pion/pionDetectorSystem.gdml")
    tar.add(sourceDir+"/geometry/pion/Lucite/pionDetectorLucite_world.gdml" ,arcname="geometry/pion/Lucite/pionDetectorLucite_world.gdml")
    tar.add(sourceDir+"/geometry/pion/Lucite/pionLuciteMaterials.xml" ,arcname="geometry/pion/Lucite/pionLuciteMaterials.xml")
    tar.add(sourceDir+"/geometry/pion/Lucite/pionLuciteMatrices.xml" ,arcname="geometry/pion/Lucite/pionLuciteMatrices.xml")
    tar.add(sourceDir+"/geometry/pion/Lucite/pionDetectorLucite.gdml" ,arcname="geometry/pion/Lucite/pionDetectorLucite.gdml")
    tar.add(sourceDir+"/geometry/pion/Lucite/pionDetectorLuciteSector.gdml" ,arcname="geometry/pion/Lucite/pionDetectorLuciteSector.gdml")
    tar.add(sourceDir+"/geometry/pion/GEM/pionDetectorGEM.gdml" ,arcname="geometry/pion/GEM/pionDetectorGEM.gdml")
    tar.add(sourceDir+"/geometry/pion/GEM/pionDetectorGEMOpenSector.gdml" ,arcname="geometry/pion/GEM/pionDetectorGEMOpenSector.gdml")
    tar.add(sourceDir+"/geometry/pion/TS/pionDetectorTS.gdml" ,arcname="geometry/pion/TS/pionDetectorTS.gdml")
    tar.add(sourceDir+"/geometry/pion/TS/pionDetectorTSOpenSector.gdml" ,arcname="geometry/pion/TS/pionDetectorTSOpenSector.gdml")
    tar.add(sourceDir+"/geometry/beampipe/downstream/beampipeDSMother.gdml" ,arcname="geometry/beampipe/downstream/beampipeDSMother.gdml")
    tar.add(sourceDir+"/geometry/beampipe/premoller/beampipeRaster.gdml" ,arcname="geometry/beampipe/premoller/beampipeRaster.gdml")
    tar.close()
      
if __name__ == '__main__':
    main()
