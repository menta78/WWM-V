# author: Lorenzo Mentaschi
import os
import unittest
import netCDF4 as nc
import numpy as np
from datetime import datetime

from wwmTestUtil import wwmTestUtil as util


class test2DExplicitPropagation(util.wwmTestTemplate):

  def test(self):
    # first run: cold start, writing output after 12 h and after 24 h
    os.system('ln -sf param_cold.nml param.nml; ln -sf wwminput_cold.nml wwminput.nml')
    launchCmd = util.getScWwmLaunchCommand(nProcesses=2)
    exitCode = os.system(launchCmd)
    if exitCode != 0:
      self.fail('schismWWM did not end correctly. Failing')
    # combining the output
    combineCmd = util.getCombineCommand(bgnParam=2, endParam=2)
    exitCode = os.system(combineCmd)
    if exitCode != 0:
      self.fail('combine_outputXX did not end correctly. Failing')
    # combining the hotstart file
    combineCmd = util.getCombineHotstartCommand(iTimeStep=108)
    exitCode = os.system(combineCmd)
    if exitCode != 0:
      self.fail('combine_hotstartXX did not end correctly. Failing')

    # moving the output of the test run
    os.system('cp -r outputs outputs_cold')

    ### CHECKING THE 1st HOTSTART FILE #
    dshf = nc.Dataset('wwm_hot.nc')
    # this simulation was run with LCYCLEHOT==F, so it should contain all the entries
    tms = nc.num2date(dshf.variables['ocean_time'][:], dshf.variables['ocean_time'].units)
    self.assertEqual(5, len(tms))
    self.assertEqual(datetime(2000, 1, 1), tms[0])
    self.assertEqual(datetime(2000, 1, 1, 6, 0, 0), tms[1])
    self.assertEqual(datetime(2000, 1, 2), tms[-1])
    self.assertEqual((5, 352, 8, 5), dshf.variables['ac'].shape)
    dshf.close()



    # second run: running from hour 12 with restart.
   #os.system('rm fort.* *.bin *.dat *.site wwmcheck.nml')
    os.system('rm outputs/schout_*.nc')
    os.system('ln -sf param_hot.nml param.nml; ln -sf wwminput_hot.nml wwminput.nml')
    os.system('ln -sf outputs_cold/hotstart_it=108.nc hotstart.nc')
    launchCmd = util.getScWwmLaunchCommand(nProcesses=2)
    exitCode = os.system(launchCmd)
    if exitCode != 0:
      self.fail('schismWWM did not end correctly. Failing')
    # combining the output
    combineCmd = util.getCombineCommand(bgnParam=2, endParam=2)
    exitCode = os.system(combineCmd)
    if exitCode != 0:
      self.fail('combine_outputXX did not end correctly. Failing')

    ## CHECKING THE OUTPUTS of the 2 runs: should be identical ###
    ##   only time and Hs should suffice  ########################
    ds0 = nc.Dataset('outputs_cold/schout_2.nc')
    hsNoHot = ds0.variables['WWM_1'][-1,:]
    tmNoHot = ds0.variables['time'][:]
    ds0.close()

    ds1 = nc.Dataset('outputs/schout_2.nc')
    hsHot = ds1.variables['WWM_1'][-1,:]
    tmHot = ds1.variables['time'][:]
    ds1.close()

    np.testing.assert_almost_equal(hsNoHot, hsHot)
    np.testing.assert_almost_equal(tmNoHot, tmHot)

    ### CHECKING THE 2nd HOTSTART FILE #
    dshf = nc.Dataset('wwm_hot.nc')
    # this simulation was run with LCYCLEHOT==T, so it should contain only the last 2 entries
    tms = nc.num2date(dshf.variables['ocean_time'][:], dshf.variables['ocean_time'].units)
    self.assertEqual(2, len(tms))
    self.assertEqual(datetime(2000, 1, 2, 12, 0, 0), tms[0])
    self.assertEqual(datetime(2000, 1, 2, 6, 0, 0), tms[1])
    self.assertEqual((2, 352, 8, 5), dshf.variables['ac'].shape)
    dshf.close()
    


    


if __name__ == '__main__':
  unittest.main()
