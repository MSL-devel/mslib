import os
import subprocess
import datetime
import time


buildTargets = ['testAtomGroup', 'testBBQ2', 'testCCD', 'testTree']


########################################################
# This function will attempt to build all of the targets
# supplied in the targets list.  It will return the
# name of the first target that fails to build.  If all
# build succesfully, then it will return an empty string.
def attemptToBuildTargets(dir, targets, numProcesses=1):
    os.chdir(dir)
    for currTarget in targets:
        command = 'make -j ' + str(numProcesses) + ' bin/' + currTarget
        retCode = subprocess.call(command, shell=True)
        if(retCode != 0):
            return currTarget

    return ''

#######################################################
# This function will edit the Makefile in the given
# dir to work with 64 bit machines.
def editMakefileFor64bit(mslDir):
    origMakeFileName = os.path.join(mslDir,'Makefile')
    newMakeFileName = os.path.join(mslDir,'Makefile.Orig')
    os.rename(origMakeFileName, newMakeFileName)

    makeFile = open(origMakeFileName, 'w')

    for line in open(newMakeFileName):
        if(line.find('#32BIT=F') == 0):
            line = line[1:]
        elif(line.find('EXTERNAL_LIB_DIR=/usr/lib') == 0):
            line = '#' + line
        elif(line.find('#EXTERNAL_LIB_DIR=/library/sharedlibs64') == 0):
            line = line[1:]
        makeFile.write(line)

    makeFile.close()



def editReleaseFile(mslDir):

    origReleaseFileName = os.path.join(mslDir,'src/release.h')
    newReleaseFileName = os.path.join(mslDir,'src/release.h.Orig')
    os.rename(origReleaseFileName, newReleaseFileName)

    relFile = open(origReleaseFileName, 'w')

    for line in open(newReleaseFileName):
        if(line.find('#define MSLVERSION') == 0):
            line = '#define MSLVERSION \"auto-'+datetime.datetime.now().strftime("%Y%m%dT%H%M%S")+'\"'
        elif(line.find('#define MSLDATE') == 0):
            line = '#define MSLDATE \"'+datetime.datetime.now().strftime("%Y%m%d")+'\"'

        relFile.write(line)

    relFile.close()
