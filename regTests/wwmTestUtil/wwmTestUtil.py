# WWM regression tests
# author: Lorenzo Mentaschi

import unittest, os, sys

mpirunEnvVar = 'SCWWM_MPIRUN'
scWwmExecutableEnvVar = 'SCWWM_EXE'
scWwmCombineOutputEnvVar = 'SCWWM_COMBOUTPUT_EXE'
scWwmCombineHotstartEnvVar = 'SCWWM_COMBHOTSTART_EXE'


def getScWwmLaunchCommand(nProcesses=2):
  mpiexe = getMpirunExecutable()
  schismWwmExe = getSchismWWMExecutable()
  launchCmd = 'mkdir outputs; ' + mpiexe + ' -np ' + str(nProcesses) + ' ' + schismWwmExe
  print('  launch command: \n' + launchCmd)
  return launchCmd


def getCombineCommand(bgnParam=1, endParam=2):
  combineExe = getCombineOutputExecutable()
  combineCmd = 'cd outputs; ' + combineExe + ' -b ' + str(bgnParam) + ' -e ' + str(endParam)
  print('  output combine command: \n' + combineCmd)
  return combineCmd


def getCombineHotstartCommand(iTimeStep=-1):
  combineExe = getCombineHotfileExecutable()
  combineCmd = 'cd outputs; ' + combineExe + ' -i ' + str(iTimeStep)
  print('  hotstart combine command: \n' + combineCmd)
  return combineCmd


def getMpirunExecutable():
  # returns the path of the mpirun executable
  if not mpirunEnvVar in os.environ:
    raise Exception(__getEnvironVarErrorMessage())
  return os.environ[mpirunEnvVar]


def getSchismWWMExecutable():
  # returns the path of the schismWWM executable
  if not scWwmExecutableEnvVar in os.environ:
    raise Exception(__getEnvironVarErrorMessage())
  return os.environ[scWwmExecutableEnvVar]


def getCombineOutputExecutable():
  # returns the path of the combineoutput executable
  if not scWwmCombineOutputEnvVar in os.environ:
    raise Exception(__getEnvironVarErrorMessage())
  return os.environ[scWwmCombineOutputEnvVar]


def getCombineHotfileExecutable():
  # returns the path of the combinehotfile executable
  if not scWwmCombineHotstartEnvVar in os.environ:
    raise Exception(__getEnvironVarErrorMessage())
  return os.environ[scWwmCombineHotstartEnvVar]
  


class wwmTestTemplate(unittest.TestCase):
  # base class for the automated regression tests

  def getTestPath(self):
    # returns the path of the test currently running
    mdl = sys.modules[self.__class__.__module__]
    path = os.path.abspath(mdl.__file__)
    dirName = os.path.dirname(path)
    return dirName

  def getRunDir(self):
    # returns the run directory of the current test
    # this is the directory where schism/wwm should run
    return os.path.join(self.getTestPath(), 'run')

  def prepareRunDirectory(self):
    # creating the directory "run" and linking all the needed files
    tstpth = self.getTestPath()
    confDir = os.path.join(tstpth, 'conf')
    runDir = self.getRunDir()
    cmd = 'rm -rf {runDir}; mkdir {runDir}; cd {runDir}; ln -s {confDir}/* ./'\
            .format(confDir=confDir, runDir=runDir)
    os.system(cmd)
    
  def cleanupRunDirectory(self):
    # removes the directory run
    os.system('rm -rf ' + self.getRunDir())
  
  def setUp(self):
    # this method is invoked before the run of each test
    # it prepares the run directory on the basis of the conf directory
    print('  running the tests in ' + self.__class__.__module__)
    self.prepareRunDirectory()
    self.launchDir = os.getcwd()
    # setting the working directory
    os.chdir(self.getRunDir())

  def tearDown(self):
    # this method is invoked after the run of each test
    os.chdir(self.launchDir)
    try:
      # getting the test outcome: maybe there were errors of some kind
      rslt = self.defaultTestResult()
      outcome = self._outcome
      self._feedErrorsToResult(rslt, outcome.errors)
      doCleanUp = len(rslt.failures) == 0 and len(rslt.errors) == 0
    except:
      doCleanUp = True
    if doCleanUp:
      # cleaning up the run directory only if the run succeeds
      self.cleanupRunDirectory()
    else:
      print('   ### TEST FAILED. THE RUN DIRECTORY WON''T BE CLEANED UP')



def __getEnvironVarErrorMessage():
  msg ="""
    You need to specify the executables paths by setting the corresponding envirnoment variables:
    - {mpirunEnvVar}: executable of mpirun, you can set this with the command
        $ export {mpirunEnvVar}=mpirun
    - {scWwmExecutableEnvVar}: executable of schismWWM, you can set this with the command
        $ export {scWwmExecutableEnvVar}=/mypath/schismWwm.exe
    - {scWwmCombineOutputEnvVar}: executable for combining the output of schism. 
        you can set this with the command
        $ export {scWwmCombineOutputEnvVar}=/mypath/combine_outputXX
    - {scWwmCombineHotstartEnvVar}: executable for combining the hotstart files of schism, needed by some tests. 
        you can set this with the command
        $ export {scWwmCombineHotstartEnvVar}=/mypath/combine_hotstartXX
    You may also need to set properly the LD_LIBRARY_PATH variable, with mpi and netcdf stuff.

    Suggestion: exports all the variables needed in the file set_envoron.sh,
    which is sourced before launching the test.
    """.format(mpirunEnvVar=mpirunEnvVar,
               scWwmExecutableEnvVar=scWwmExecutableEnvVar, 
               scWwmCombineOutputEnvVar=scWwmCombineOutputEnvVar,
               scWwmCombineHotstartEnvVar=scWwmCombineHotstartEnvVar)
  return msg

