import mslBuildTools
import sys
import miscUtils
import getpass
import os

OPTIONS = {'testlevel' :                  (True, 'The test level file.'),
           'numprocesses' :               (False, 'How many processes should we use while compiling?', 1),}

SHARED_SUBMIT_DIR = '/tmp'

################################################################################
args = miscUtils.check_options(sys.argv, OPTIONS)

#currUser = getpass.getuser()
#newDirName = miscUtils.getRandomFileName(os.path.join(SHARED_SUBMIT_DIR, currUser))
#newMslDirName = os.path.join(newDirName, 'trunk')
newMslDirName = './'

results = mslBuildTools.run_tests(newMslDirName, mslBuildTools.buildTargets, int(args['numprocesses']))
mslBuildTools.print_test_results(results)
