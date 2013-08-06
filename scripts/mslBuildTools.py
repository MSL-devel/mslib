import os
import subprocess
import datetime
import time


buildTargets = 'testLevels/submit.level'

RESULTS_FILE = 'RESULTS.txt'

########################################################
# This function will attempt to build all of the targets
# supplied in the targets list.  It will return the
# name of the first target that fails to build.  If all
# build succesfully, then it will return an empty string.
def attempt_to_build_targets(dir, targets, numProcesses=1,boost=False, gsl=False, glpk=False, R=False):
    os.chdir(dir)

    #export_command = 'export MSL_R=F;export MSL_BOOST=F;export MSL_GSL=F; export MSL_GLPK=F;'
    export_command = ''
    libs_on="lib"
    if boost:
        export_command += 'export MSL_BOOST=T;'
        libs_on += "B"
    else:
        export_command += 'export MSL_BOOST=F;'
    if gsl:
        export_command += 'export MSL_GSL=T;'
        libs_on += "G"
    else:
        export_command += 'export MSL_GSL=F;'
    if glpk:
        export_command += 'export MSL_GLPK=T;'        
        libs_on += "K"
    else:
        export_command += 'export MSL_GLPK=F;'        
    if R:
        export_command += 'export MSL_R=T;'        
        libs_on += "R"
    else:
        export_command += 'export MSL_R=F;'        
    
    try:
        command = export_command + 'make clean > /dev/null;'
        retCode = subprocess.call(command, shell=True)
        
        for currTarget in targets:
            command = export_command + 'make -j ' + str(numProcesses) + ' ' + currTarget +  ('  2>> compile.%s.err | tee -a compile.%s.out' % (libs_on,libs_on))
            retCode = subprocess.call(command, shell=True)
            if (retCode != 0):
                return currTarget
    except:
        print 'Exception in building targets.'

    return ''

########################################################
# This function will read in tests from a test file and
# return a dictionary with the target as the key
# and the command line as the value.
def read_test_file(testFile, seenFiles = {}):
    tests = {}
    seenFiles[testFile] = 1
    
    for l in open(testFile):
        l = l.strip()
        if((l != '') and (l[0] != '#')):
            w = l.split()
            # We want to prevent hanging on circular dependencies,
            # so if we encounter an import statement for a file
            # we've already supposedly imported, skip over it.
            if( (w[0] == 'import') and (w[1] not in seenFiles)):
                newTests = read_test_file(w[1], seenFiles)
                tests.update(newTests)
            else:
                tests[w[0]] = (' '.join(w[1:-1]), w[-1])
    
    return tests

########################################################
# This function read through the text output of a test
# and look for the words LEAD or GOLD on a single line.
# It will return the last instance of either word.
def check_for_test_status(resultFile):
    status = ''
    
    for l in open(resultFile):
        l = l.strip()
        if( (l == 'GOLD') or (l == 'LEAD') or (l == 'PASS') or (l == 'FAIL')):
            status = l
    
    return status

########################################################
# This function read through the text output of a test
# and look for the words LEAD or GOLD on a single line.
# It will return the last instance of either word.
def print_test_results(results, skipFailures = True):
    print 'TEST RESULTS:'
    for result in results:
        if(not skipFailures or (result != 'failures')):
            print '  ' + result + '\t' + results[result][0]

########################################################
# This funciton will just check to see that each test
# passed the minimum required for it.
def check_results_for_failures(testResults):
    for test in testResults:
        if(test != 'failures'):
            
            # A test passes if it is the same value as the expected value,
            # if it is GOLD, or if it is LEAD and there was no expected value.
            if( (testResults[test][0] == testResults[test][1]) or 
                (testResults[test][0] == 'GOLD') or 
                ((testResults[test][0] == 'LEAD') and (testResults[test][1] == '')) ):
                continue
            else:
                testResults['failures'] += [test]

########################################################
# This function will read in tests from a test file,
# build them, run them, and return a dictionary indicating
# the test result.
# AS 20130414 Turned GSL on by default
def run_tests(dir, testFile, numProcesses=1,boost=False, gsl=True, glpk=False, R=False):
    testResults = {}
    testResults['failures'] = []
    
    tests = read_test_file(testFile)
    
    failures = attempt_to_build_targets(dir, tests.keys(), numProcesses, boost, gsl, glpk, R)
    if(failures != ''):
        testResults['failures'] += [failures]
        return testResults
    
    os.chdir(dir)
    for test in tests:
        if(os.path.exists(test) == False):
            testResults['failures'] += [test]

        cmdLine = tests[test][0] + ' > ' + RESULTS_FILE
        print 'Running: ' + cmdLine
        subprocess.call(cmdLine, shell=True)
        result = check_for_test_status(RESULTS_FILE)
        testResults[test] = (result, tests[test][1])
    
    check_results_for_failures(testResults)

    return testResults

#######################################################
# This function will edit the Makefile in the given
# dir to work with 64 bit machines.
def edit_makefile_for_64bit(mslDir):
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



def edit_release_file(mslDir):

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
