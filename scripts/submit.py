#!/usr/bin/python

import sys
import os
import subprocess
import getpass
import shutil
import string
import re
from optparse import OptionParser
import mslBuildTools
import miscUtils


FILE_LIST_FILE_NAME = 'svn_submit_file.txt'
SHARED_SUBMIT_DIR = '/tmp'
FAILED_SUBMIT_DIR = os.path.join(SHARED_SUBMIT_DIR, 'failed')
RELEASE_FILE = os.path.join('src','release.h')
SVN_ROOT_PATH = 'svn+ssh://%s@svn.code.sf.net/p/mslib/code/'


def get_options():
    '''This function will create the option parser for the submit script.'''
    parser = OptionParser()
    parser.add_option('-f', '--file', dest='fileList', 
                      help='Comma delimited list of files to submit.', metavar='FILES')
    try:
        parser.add_option('-u', '--user', dest='userName', default=getpass.getuser(), help='Sourceforge username.')
    except:
        parser.add_option('-u', '--user', dest='userName', default='no_name', help='Sourceforge username.')
    parser.add_option('-m', '--message', dest='fileMessage', 
                      help='Change description for the given files.', metavar='MESSAGE')
    parser.add_option('-d', '--delete', dest='filesToDelete', 
                      help='Comma delimited list of files to delete.', metavar='FILES')
    parser.add_option('-r', '--restart', dest='restart', default=False, action='store_true',
                      help='This will allow you to restart your submit filelist.')
    parser.add_option('-k', '--nukefaileddir', dest='nuke', default=False, action='store_true',
                      help='This will delete all files from your failed submit directory.')
    parser.add_option('-c', '--numcores', dest='numCores', default=1, type='int',
                      help='How many cores should we use during our build?')
    parser.add_option('-s', '--skip_build', dest='skipBuild', default=False, action='store_true',
                      help='THIS IS VERY DANGEROUS!!!!\nThis flag will skip building the tree and running tests.\nOnly use if you absolutely know what you are doing.')
    parser.add_option('-n', '--now', dest='doSubmit', default=False,
                      action='store_true', help='When you have entered all of the files you wish to add and are ready to do the submit.')
                      
    return parser

########################################################
# This function will append the filenames and messages
# given on the command line to the file used to
# actually submit.
def addFilesAndMessagesToFileListFile(fullFileListFileName, filenames, f_or_d, messages):
    fileList = ''
    for filename in filenames:
        fileList += '"' + filename + '",'
    fileList = fileList[0:-1]


    fileListFile = open(fullFileListFileName, 'a')
    fileListFile.write('["' + f_or_d + '",[' + fileList + '],"' + messages + '"]\n')
    fileListFile.close()

########################################################
# This function will create and sync the given dir.
def createAndSyncDir(newDirName, mslSubDir, userName):
    os.mkdir(newDirName)
    os.chdir(newDirName)
    #command = 'cvs co msl'
    if(mslSubDir != 'trunk'):
        mslSubDir = 'branches/' + mslSubDir
    command = 'svn co ' + SVN_ROOT_PATH + '/' + mslSubDir
    subprocess.call(command % (userName),shell=True)

########################################################
# This function will open the given file name,
# and take note of all of the files that will be modified.
def getFilesToModify(fullFileListFileName):
    filesToModify = []
    
    for line in open(fullFileListFileName):
        if(line.strip() != ''):
            temp = eval(line.strip())
            if(temp[0] == 'f'):
                filesToModify += temp[1]

    return filesToModify

########################################################
# This function will open the given file name,
# and take note of all of the files that will be deleted.
def getFilesToDelete(fullFileListFileName):
    filesToDelete = []

    for line in open(fullFileListFileName):
        if(line.strip() != ''):
            temp = eval(line.strip())
            if(temp[0] == 'd'):
                filesToDelete += temp[1]

    return filesToDelete

########################################################
# This function will submit the files.
def submitFiles(cwd, newMslDirName, dirsToBeAdded, fileListFileName, newVersion, releaseFileName, mslSubDir, userName):
    os.chdir(newMslDirName)
    for dirToAdd in dirsToBeAdded:
        command = 'svn commit -N -m "Adding directory." ' + dirToAdd
        subprocess.call(command, shell=True)
        
    for line in open(fileListFileName):
        temp = eval(line.strip())
        command = 'svn commit -m "' + temp[2]

        if(temp[0] == 'f'):
            command += '"'
            for currFile in temp[1]:
                command += ' ' + currFile
        elif(temp[0] == 'd'):
            for currFile in temp[1]:
                command += " '" + currFile + "'"
            command += '"'

        print 'Commiting with the following command: ' + command
        subprocess.call(command, shell=True)

    command = 'svn commit -m "Rolling release file." ' + releaseFileName
    print 'Commiting with the following command: ' + command
    subprocess.call(command, shell=True)

    if(mslSubDir != 'trunk'):
        mslSubDir = 'branches/' + mslSubDir
    command = 'svn copy -m "Copying version." ' + SVN_ROOT_PATH + '/' + mslSubDir + ' ' + SVN_ROOT_PATH + 'tags/' + newVersion
    command = command % (userName, userName)
    print 'Copying repository with the following command: ' + command
    subprocess.call(command, shell=True)
    
    os.chdir(cwd)


########################################################
# Copy the files from the current tree to the test tree.
def copyMyFiles(cwd, mslDir, myFiles):
    for currFile in myFiles:
        srcFullFilePath = os.path.join(cwd, currFile)
        dstFullFilePath = os.path.join(mslDir, currFile)
        shutil.copy2(srcFullFilePath, dstFullFilePath)

########################################################
# This function will look to see what files and dirs
# are new, and therefore need to be added.
def addFilesAndDirectories(cwd, newMslDirName, myFiles):
    errFileName = 'submit.err'
    dirsToBeAdded = []
    filesToBeAdded = []
    os.chdir(newMslDirName)

    for currFile in myFiles:
        # First test to see if the directory even exists.
        # If not, add it.
        dirName = os.path.split(currFile)[0]
        fullDirName = os.path.join(newMslDirName,dirName)

        if(os.path.isdir(fullDirName) == False):
            os.mkdir(fullDirName)
            #command = 'cvs add ' + dirName
            dirsToBeAdded += [dirName]
            command = 'svn add ' + dirName
            subprocess.call(command, shell=True)

        # Now check to see if the file exists.  
        if(os.path.isfile(currFile) == False):
            filesToBeAdded += [currFile]

    # Change back to the original directory
    os.chdir(cwd)
    return (dirsToBeAdded, filesToBeAdded)

########################################################
# This function will look in the release file to 
# find the release version number.
def getVersionNumber(releaseFilename):
    version = ''
    for line in open(releaseFilename):
        line = line.strip()
        if(line.find('MSLVERSION') > 0):
            words = line.split()
            if(len(words) >= 3):
                version = words[2].strip('"')

    return version

    
########################################################
# This function will add the given files.
def addFiles(cwd, newMslDirName, filesToBeAdded):
    os.chdir(newMslDirName)

    for currFile in filesToBeAdded:
        #command = 'cvs add ' + currFile
        command = 'svn add ' + currFile
        subprocess.call(command, shell=True)

    os.chdir(cwd)

########################################################
# This function will delete files and directories.
def deleteFiles(cwd, newMslDirName, filesToBeDeleted):
    os.chdir(newMslDirName)
    
    for currFile in filesToBeDeleted:
        command = 'svn delete ' + currFile
        subprocess.call(command, shell=True)
    
    os.chdir(cwd)

########################################################
# This function will simply take a long line of text
# and split it into smaller lines for the release file.
def createTextForReleaseFile(line):
    out = ''
    
    nl = '               '
    for word in line.split():
        nl += ' ' + word
        if(len(nl) > 125):
            out += nl + '\n'
            nl = '                '

    return out + nl

########################################################
# This function will modify the release.h file.
def modifyReleaseFile(releaseFileName, newVersionNumber, userName, fileListFileName, skipBuild):
    import datetime
    import shutil
    now = datetime.date.today().strftime('%B %d, %Y')
    
    tempFileName = miscUtils.getRandomFileName(releaseFileName)
    tempFile = open(tempFileName, 'w')
    
    seenHistory = False
    for line in open(releaseFileName):
        if(seenHistory == False):
            if(line.find('MSLVERSION') >= 0):
                line = '#define MSLVERSION "' + newVersionNumber + '"\n'
            if(line.find('MSLDATE') >= 0):
                line = '#define MSLDATE "' + now + '"\n'
        
        tempFile.write(line)
        
        if(line.find('HISTORY:') >= 0):
            templine = newVersionNumber + '    ' + now + '    ' + userName + '\n'
            if(skipBuild):
                templine += '               WARNING!  Files submitted without building tree or running tests.\n'
            
            for submitline in open(fileListFileName):
                #temp = eval(submitline.strip())[1]
                temp = eval(submitline.strip())
                files = str(temp[1]).strip('[]')
                while(len(files) < 15):
                    files = files + ' '
                comment = temp[2]
                temp = createTextForReleaseFile(files + ' -' + comment) + '\n'
                templine += temp
            tempFile.write(templine)
            seenHistory = True
            
    tempFile.close()
    
    shutil.copyfile(tempFileName, releaseFileName)
    os.remove(tempFileName)

########################################################
# This function will just make sure that the new
# version is greater than the old version.
def isNewVersionLegal(newVersion, myVersion):
    newVersion = newVersion.split('.')
    myVersion = myVersion.split('.')

    nv = [int(b) for b in newVersion]
    mv = [int(b) for b in myVersion]

    if(len(nv) != len(mv)):
        return False

    for i in range(len(mv)):
        if (nv[i] == mv[i]):
            continue
        if(nv[i] < mv[i]):
            return False
        if(nv[i] > mv[i]):
            return True

    return False

########################################################
# Main
parser = get_options()
(options, args) = parser.parse_args()

cwd = os.getcwd()
# Append the files to the filelist
fullFileListFileName = os.path.join(cwd, FILE_LIST_FILE_NAME)
mslSubDir = os.path.split(cwd)[1]

if(options.restart):
    if(os.path.isfile(fullFileListFileName)):
        os.remove(fullFileListFileName)

if(options.nuke):
    command = 'rm -rf ' + FAILED_SUBMIT_DIR
    subprocess.call(command,shell=True)


if((options.doSubmit == False)):
    print 'Note, submit is just queing the files for submit.'
    print 'To actually continue the submit, please specify the --now option'

if(options.fileList != None):
    options.fileList = re.split('[, ]', options.fileList)
    addFilesAndMessagesToFileListFile(fullFileListFileName, options.fileList, 'f', options.fileMessage)

if(options.filesToDelete != None):
    options.filesToDelete = re.split('[, ', options.filesToDelete)
    addFilesAndMessagesToFileListFile(fullFileListFileName, options.filesToDelete, 'd', options.fileMessage)

# If now was specified on the command line,
# then we want to actually do the submit now.
releaseFileModified = False
try:
    if(options.doSubmit):
        myFiles = getFilesToModify(fullFileListFileName)

        newDirName = miscUtils.getRandomFileName(os.path.join(SHARED_SUBMIT_DIR, options.userName))
        createAndSyncDir(newDirName, mslSubDir, options.userName)
        #newMslDirName = os.path.join(newDirName,'msl')
        newMslDirName = os.path.join(newDirName, mslSubDir)
        # Set the MSL_DIR environmental directory to the new, test trunk.
        os.environ['MSL_DIR'] = newMslDirName

        myVersion = getVersionNumber(os.path.join(cwd, RELEASE_FILE))
        topVersion = getVersionNumber(os.path.join(newMslDirName, RELEASE_FILE))
        if(myVersion != topVersion):
            print 'The version given in your release.h, ' + myVersion + ', does not match that at the top of the tree, ' + topVersion + '.'
            print 'Please sync your tree to the top of the tree, merge any conflicting files, and try to submit again.'
            print 'Your list of files to submit will be kept, so once the tree is synced and merged, just run submit.py -now again.'
            sys.exit(1)

        versionLegal = False
        while(versionLegal == False):
            newVersion = raw_input('The current version number is: ' + myVersion + '\nWhat is the new version number?\n')
            versionLegal = isNewVersionLegal(newVersion, myVersion)
            if(versionLegal == False):
                print 'The new version of ' + newVersion + ' is not legal.  Try again.\n'

        # Modify the release.h file and and add it to the list.
        modifyReleaseFile(os.path.join(cwd, RELEASE_FILE), newVersion, options.userName, fullFileListFileName, options.skipBuild)
        releaseFileModified = True
        myFiles += [RELEASE_FILE]
        filesToBeDeleted = getFilesToDelete(fullFileListFileName)

        shutil.move(fullFileListFileName, newDirName)
        mslBuildTools.edit_makefile_for_64bit(newMslDirName)

        (dirsToBeAdded, filesToBeAdded) = addFilesAndDirectories(cwd, newMslDirName, myFiles)
        copyMyFiles(cwd, newMslDirName, myFiles)
        addFiles(cwd, newMslDirName, filesToBeAdded)
        deleteFiles(cwd, newMslDirName, filesToBeDeleted)
        
        results = {'failures': []}
        if(options.skipBuild == False):
            results = mslBuildTools.run_tests(newMslDirName, mslBuildTools.buildTargets, options.numCores)
        
        # Currently aren't checking if tests are passing, failing, lead, or gold.
        if(len(results['failures']) == 0):
            mslBuildTools.print_test_results(results)
            submitFiles(cwd, newMslDirName, dirsToBeAdded, os.path.join(newDirName, FILE_LIST_FILE_NAME), newVersion, RELEASE_FILE, mslSubDir, options.userName)
            subprocess.call('rm -rf ' + newDirName, shell=True)
            print 'Submitted!'
        else:
            shutil.move(newDirName, FAILED_SUBMIT_DIR)
            newDir = os.path.join(FAILED_SUBMIT_DIR, os.path.split(newDirName)[1])
            print 'Submit failed on target: ' + ' '.join(results['failures']) + '.  Moved directory to ' + newDir + '.'
            if(releaseFileModified):
                print "A test program failed"	
                subprocess.call('svn revert ' + os.path.join(cwd, RELEASE_FILE), shell=True)
            sys.exit(1)
except Exception, e:
    if(releaseFileModified):
        print "An error has occurred: %s" % e
        subprocess.call('svn revert ' + os.path.join(cwd, RELEASE_FILE), shell=True)
