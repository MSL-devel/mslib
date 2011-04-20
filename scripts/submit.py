import miscUtils
import sys
import os
import subprocess
import getpass
import shutil
import mslBuildTools
import string

fileListFileName = 'svn_submit_file.txt'
sharedSubmitDir = '/tmp'
failedSubmitDir = os.path.join(sharedSubmitDir, 'failed')
releaseFile = os.path.join('src','release.h')


########################################################
# This will check to see that the proper arguments
# were given and print a usage message if not.
def testOptionsAndPrintUsage(options):
    filesAndMessageSet = ('f' in options) and ('m' in options)
    nowSetAndFileExists = ('now' in options) and os.path.isfile(fileListFileName)
    restartOrNukeExists = ('restart' in options) or ('nukefaileddir' in options)

    if( (filesAndMessageSet or nowSetAndFileExists or restartOrNukeExists) == False):
        print 'Usage: submit -f <list of files> -m <file descriptions>'
        print '       submit -now (when you have entered all of the files you wish to add.)'
        print '   or  submit -restart (this will allow you to restart your submit filelist)'
        print '   or  submit -nukefaileddir (this will delete all files from your failed directory'
        sys.exit(1)

########################################################
# This function will append the filenames and messages
# given on the command line to the file used to
# actually submit.
def addFilesAndMessagesToFileListFile(fullFileListFileName, filenames, messages):
    fileList = ''
    for filename in filenames:
        fileList += '"' + filename + '",'
    fileList = fileList[0:-1]


    fileListFile = open(fullFileListFileName, 'a')
    fileListFile.write('[[' + fileList + '],"' + string.join(messages) + '"]\n')
    fileListFile.close()

########################################################
# This function will create and sync the given dir.
def createAndSyncDir(newDirName):
    os.mkdir(newDirName)
    os.chdir(newDirName)
    #command = 'cvs co msl'
    command = 'svn co https://mslib.svn.sourceforge.net/svnroot/mslib/trunk'
    subprocess.call(command,shell=True)

########################################################
# This function will open the given file name,
# and take note of all of the files that will be modified.
def getFilesToModify(fullFileListFileName):
    filesToModify = []
    
    for line in open(fullFileListFileName):
        if(line.strip() != ''):
            temp = eval(line.strip())
            filesToModify += temp[0]

    return filesToModify

########################################################
# This function will submit the files.
def submitFiles(newMslDirName, fileListFileName, newVersion, releaseFileName):
    for line in open(fileListFileName):
        temp = eval(line.strip())
        #command = 'cvs commit -m "' + temp[1] + '" '
        command = 'svn commit -m "' + temp[1] + '"'

        for currFile in temp[0]:
            command += ' ' + currFile

        print 'Commiting with the following command: ' + command
        subprocess.call(command, shell=True)

    command = 'svn commit -m "Rolling release file." ' + releaseFileName
    print 'Commiting with the following command: ' + command
    subprocess.call(command, shell=True)

    command = 'svn copy -m "Copying version." https://mslib.svn.sourceforge.net/svnroot/mslib/trunk https://mslib.svn.sourceforge.net/svnroot/mslib/tags/' + newVersion
    print 'Copying repository with the following command: ' + command
    subprocess.call(command, shell=True)


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
            command = 'svn add ' + dirName
            subprocess.call(command, shell=True)

        # Now check to see if the file exists.  
        if(os.path.isfile(currFile) == False):
            filesToBeAdded += [currFile]

    # Change back to the original directory
    os.chdir(cwd)
    return filesToBeAdded

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
def modifyReleaseFile(releaseFileName, newVersionNumber, userName, fileListFileName):
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
            
            for submitline in open(fileListFileName):
                #temp = eval(submitline.strip())[1]
                temp = eval(submitline.strip())
                files = str(temp[0]).strip('[]')
                while(len(files) < 15):
                    files = files + ' '
                comment = temp[1]
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
options = miscUtils.parseCommandLineOptions(sys.argv[1:])
cwd = os.getcwd()
# Append the files to the filelist
fullFileListFileName = os.path.join(cwd, fileListFileName)

if('restart' in options):
    if(os.path.isfile(fullFileListFileName)):
        os.remove(fullFileListFileName)

if('nukefaileddir' in options):
    command = 'rm -rf ' + failedSubmitDir
    subprocess.call(command,shell=True)

testOptionsAndPrintUsage(options)


if(('now' not in options) and ('f' in options)):
    print 'Note, submit is just queing the files for submit.'
    print 'To actually continue the submit, please specify the -now option'

if('f' in options):
    addFilesAndMessagesToFileListFile(fullFileListFileName, options['f'], options['m'])

# If now was specified on the command line,
# then we want to actually do the submit now.
if('now' in options):
    myFiles = getFilesToModify(fullFileListFileName)

    currUser = getpass.getuser()
    newDirName = miscUtils.getRandomFileName(os.path.join(sharedSubmitDir,currUser))
    createAndSyncDir(newDirName)
    #newMslDirName = os.path.join(newDirName,'msl')
    newMslDirName = os.path.join(newDirName,'trunk')
    
    myVersion = getVersionNumber(os.path.join(cwd, releaseFile))
    topVersion = getVersionNumber(os.path.join(newMslDirName, releaseFile))
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
    modifyReleaseFile(os.path.join(cwd, releaseFile), newVersion, currUser, fullFileListFileName)
    myFiles += [releaseFile]

    shutil.move(fullFileListFileName, newDirName)
    mslBuildTools.editMakefileFor64bit(newMslDirName)

    filesToBeAdded = addFilesAndDirectories(cwd, newMslDirName, myFiles)
    copyMyFiles(cwd, newMslDirName, myFiles)
    addFiles(cwd, newMslDirName, filesToBeAdded)

    failedBuilds = 'submit.py'
    try:
        numProcesses = int(options['now'][0])
        failedBuilds = mslBuildTools.attemptToBuildTargets(newMslDirName, mslBuildTools.buildTargets, numProcesses)
    except:
        failedBuilds = mslBuildTools.attemptToBuildTargets(newMslDirName, mslBuildTools.buildTargets)

    if(failedBuilds == ''):
        submitFiles(newMslDirName, os.path.join(newDirName, fileListFileName), newVersion, releaseFile)
        subprocess.call('rm -rf ' + newDirName, shell=True)
        print 'Submitted!'
    else:
        shutil.move(newDirName, failedSubmitDir)
        newDir = os.path.join(failedSubmitDir,os.path.split(newDirName)[1])
        print 'Submit failed on target: ' + failedBuilds + '.  Moved directory to ' + newDir + '.'
        sys.exit(1)
