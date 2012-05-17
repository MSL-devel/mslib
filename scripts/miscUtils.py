import glob
import os
import random

################################################
# The user passes in a filename, and this
# function simply forms a list with an entry
# per line of the file.  It will ignore lines
# where the first non-whitespace char is a #.
def loadListFromFile(fileName):
    myList = []
    
    for line in open(fileName):
        line = line.strip()
        # Any line that has the first non-whitespace character
        # as a # sign is treated as a comment.
        if((len(line) > 0) and (line[0] != '#')):
            myList += [line]
    
    return myList

################################################
# The user passes in either a directory or filename.
# If a directory, then this function will return
# all files in the given directory matching the
# query sequence.  By default it will recursively
# walk through all subdirectories.  If a file,
# The function will simply extract all filenames
# listed in the given file.
def getFileList(dirOfInterestOrFile, query='*', recursive=True):
    listOfFiles = []
    # Recursively get all files in this directory and
    # sub-directory.
    if(os.path.isdir(dirOfInterestOrFile)):
        if(recursive):
            for path, dirs, files in os.walk(dirOfInterestOrFile):
                listOfFiles += glob.glob( os.path.join(path, query) )
        else:
            listOfFiles += glob.glob( os.path.join(dirOfInterestOrFile, query) )
    elif(os.path.isfile(dirOfInterestOrFile)):
        for line in open(dirOfInterestOrFile):
            listOfFiles += [line.strip()]

    return listOfFiles

################################################
# The user passes in either a directory or filename.
# If a directory, then this function will return
# all directories in the given directory.
# By default it will recursively
# walk through all subdirectories.  If a file,
# The function will simply extract all directories
# listed in the given file.
def getDirList(dirOfInterestOrFile, recursive=True):
    listOfDirs = []
    # Recursively get all directories in this directory and
    # sub-directory.
    if(os.path.isdir(dirOfInterestOrFile)):
        if(recursive):
            for path, dirs, files in os.walk(dirOfInterestOrFile):
                listOfDirs += [path]
        else:
            temp = os.listdir(dirOfInterestOrFile)
            for currTest in temp:
                if(os.path.isdir(currTest)):
                    listOfDirs += [currTest]
    elif(os.path.isfile(dirOfInterestOrFile)):
        for line in open(dirOfInterestOrFile):
            listOfDirs += [line.strip()]

    return listOfDirs


################################################
# This function will return a unique filename with
# the given prefix, followed by some random digits
def getRandomFileName(prefix):

    fileName = ''
    fileExists = True

    # Loop until we find a unique filename.
    while(fileExists):
        # Choose some random 4 byte hex number.
        randNum = hex(random.randint(0, 2**32 - 1))[2:]
        # Append the random number to the prefix.
        fileName = prefix + '_' + randNum
        # Does this file already exist?
        if(not (os.path.isfile(fileName) or os.path.isdir(fileName))):
            fileExists = False

    # Return unique filename.
    return fileName

##############################################
# This function will poll on the given filename
# looking for the given key word.  It will return
# once it has seen the word.
def pollOutputFile(outputFilename, keyWord):
    while (not os.path.exists(outputFilename)):
        a = 1

    # Open file.
    myFile = open(outputFilename)
    notDone = True

    # Loop until we see the key word.
    while ( notDone ):
        line = myFile.readline()
        if (line and line.strip() == keyWord):
            notDone = False

##############################################
# This function will poll on the given filename
# looking for the given key word.  It will return
# once it has seen the word.
def launchCommandsInParallelAndPollForDone(commandList, queue):
    magicWord = 'DONE'
    # First create file names for each of our notifier files.
    notifierFiles = []
    for index in range(len(commandList)):
        prefix = './_notifier_' + str(index)
        notifierFiles += [getRandomFileName(prefix)]

    # Now create new command list.
    newCommandList = []
    for index in range(len(commandList)):
        newCommand = queue + ' -qsub "'
        newCommand += commandList[index]
        newCommand += ' ; echo "' + magicWord + '" > ' + notifierFiles[index] +'"'
        newCommandList += [newCommand]

    # Now actually launch the commands
    for command in newCommandList:
        print 'Launching ' + command
        os.system(command)

    # Poll for done and remove notifier file
    # once we have seen it.
    for notifierFile in notifierFiles:
        print 'Polling on ' + notifierFile + '.'
        pollOutputFile(notifierFile, magicWord)
        os.remove(notifierFile)
    

##############################################
def launchCommandsInParallelAndPollForDoneLocal(commandList, numConcurrentProcess):
    import os
    import stat
    scriptFileNames = []

    for i in range(numConcurrentProcess):
        scriptFileNames += [getRandomFileName('_script_' + str(i))]
    
    currScriptIndex = 0
    while(len(commandList) > 0):
        currScriptFileName = scriptFileNames[currScriptIndex]
        if(os.path.isfile(currScriptFileName) == False):
            outFile = open(currScriptFileName, 'w')
            outFile.write('#!/bin/sh\n')
            outFile.write(commandList[0] + '\n')
            outFile.write('rm ' + currScriptFileName + '\n')
            outFile.close()
            os.chmod(currScriptFileName, stat.S_IRWXU)
            os.system(currScriptFileName + '&')
            print 'Launching ' + commandList[0] + '.'
            commandList = commandList[1:]

        currScriptIndex = (currScriptIndex+1) % len(scriptFileNames)

    filesStillExist = True
    print 'Waiting for last commands to complete'
    while(filesStillExist):
        filesStillExist = False
        for currScriptFileName in scriptFileNames:
            filesStillExist = filesStillExist or os.path.isfile(currScriptFileName)
    print 'DONE'

##############################################
# This function will loop over the arguments
# given in argv and create a dictionary with
# key value of the parameter name and value
# of the given value.
def parseCommandLineOptions(argv, convertToLower = True):
    args = {}

    # Some people use two dashes, some people use one dash.
    import string

    argv = string.join(argv)
    options = argv.strip().split('--')

    for currOption in options:
        # If the current option is empty, ignore.
        if(currOption != ''):
            words = currOption.split()
            key = words[0]
            if(convertToLower):
                key = key.lower()
            val = ['']
            if(len(words) > 1):
                val = words[1:]

            args[key] = val

    return args

###############################################
## This function will loop over the arguments
## given in argv and create a dictionary with
## key value of the parameter name and value
## of the given value.
#def parseCommandLineOptions(argv, convertToLower = True, preserveSlashDashes = True):
#    args = {}
#
#    # Some people use two dashes, some people use one dash.
#    import string
#
#    argv = string.join(argv).replace('--', '-')
#    if(preserveSlashDashes):
#        argv = argv.replace('\-', 'itsadash!')
#
#    options = argv.strip().split('-')
#
#
#    for currOption in options:
#        # If the current option is empty, ignore.
#        if(currOption != ''):
#            words = currOption.split()
#            key = words[0]
#            if(convertToLower):
#                key = key.lower()
#            val = ['']
#            if(len(words) > 1):
#                val = words[1:]
#
#            if(preserveSlashDashes):
#                temp = []
#                for w in val:
#                    w = w.replace('itsadash!', '-')
#                    temp += [w]
#                val = temp
#
#            args[key] = val
#
#    return args

##############################################
# This function will loop over the arguments
# given in argv and create a dictionary with
# key value of the parameter name and value
# of the given value.
def parseCommandLineOptionsFromFile(fn, convertToLower = True):
    args = {}
    import os
    if(os.path.exists(fn)):
        commandLineEquivalent = ''
        for l in open(fn):
            w = l.strip().split()

            if( (len(w) >= 1) and (w[0][0] != '#') ):
                commandLineEquivalent += '--' + w[0] + ' '
                if(len(w) > 1):
                    commandLineEquivalent += w[1].replace('-', '\-') + ' '

        args = parseCommandLineOptions(fn, convertToLower, True)

    return args

################################################################################
# This function will print the expected usage of the given script.
def print_usage(OPTIONS, scriptName, missingKey = ''):
    if(missingKey == ''):
        print 'Missing the required option ' + missingKey + '.'

    requiredArgs = ''
    optionalArgs = ''
    for key in OPTIONS:
        if(OPTIONS[key][0]):
            requiredArgs += '    --' + key + ' (' + OPTIONS[key][1] + ')\n'
        else:
            optionalArgs += '    --' + key + ' (' + OPTIONS[key][1] + ')\n'

    print 'Usage: python ' + scriptName
    print '  Required Arguments:'
    print requiredArgs.rstrip()
    print '  Optional Arguments:'
    print optionalArgs.rstrip()


################################################################################
# This funciton will check that we have values for all of the required options.
# It will aslo return a dictionary with the proper args.
def check_options(argv, OPTIONS, scriptName=''):
    args = {}
    import sys
    if(scriptName == ''):
        scriptName = argv[0]

    if(len(argv) == 1):
        print_usage(OPTIONS, scriptName)
        sys.exit(1)
    elif(len(argv) == 2):
        args = parseCommandLineOptionsFromFile(argv[1])
    else:
        args = parseCommandLineOptions(argv)

    for key in OPTIONS:
        if(OPTIONS[key][0] and (key not in args)):
            print 'Missing ' + key + '.'
            print_usage(OPTIONS, scriptName, key)
            sys.exit(1)
        elif( (len(OPTIONS[key]) == 3) and (key not in args)):
               args[key] = OPTIONS[key][2]
        elif(key in args):
            temp = args[key]
            # If there was no var for this arg, then it will
            # be treated as a binary variable, and set to True.
            if(temp == ['']):
                args[key] = True
            else:
                args[key] = ' '.join(temp)

    return args

##############################################
# This function will split the given list
# into numSplits different lists
def splitList(myList, numSplits):
    import math
    result = []

    numElements = len(myList)
    numElementsPerSplit = int( math.floor( float(numElements) / numSplits) )

    for i in range(numSplits-1):
        result += [myList[i*numElementsPerSplit:(i+1)*numElementsPerSplit]]

    result += [myList[(numSplits-1)*numElementsPerSplit:]]

    return result

##############################################
# This function will delete all jobs running
# on the pbs system from the current user.
def killAllPbsJobs():
    import getpass
    import subprocess

    currUser = getpass.getuser()
    tempFileName = getRandomFileName('_pbsq_info_')

    command = 'pbsq -qstat | grep ' + currUser + ' > ' + tempFileName
    subprocess.call(command, shell=True)

    pbsIds = []
    for line in open(tempFileName):
        index = line.find('.pbsserver')
        if(index > 0):
            pbsIds += [line[0:index]]

    for currId in pbsIds:
        command = 'pbsq -qdel ' + currId
        subprocess.call(command, shell=True)

    os.remove(tempFileName)

################################################
# This function will find all instances of
# the given string in a file and replace it
# with the new string.
def findAndReplaceInFile(fileName, f, r):
    import shutil
    import os

    outFileName = getRandomFileName(fileName)
    outFile = open(outFileName, 'w')
    for line in open(fileName):
        newLine = line.replace(f, r)
        outFile.write(newLine)
    outFile.close()
    shutil.copyfile(outFileName, fileName)
    os.remove(outFileName)
