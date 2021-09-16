# author: Lorenzo Mentaschi
import os
import unittest
import netCDF4 as nc
import numpy as np
from datetime import datetime

from wwmTestUtil import wwmTestUtil as util


class testUOST(util.wwmTestTemplate):

  def test(self):
    # launching schismWWM
    launchCmd = util.getScWwmLaunchCommand(nProcesses=2)
    exitCode = os.system(launchCmd)
    if exitCode != 0:
      self.fail('schismWWM did not end correctly. Failing')
    # combining the output
    combineCmd = util.getCombineCommand(bgnParam=10, endParam=10)
    exitCode = os.system(combineCmd)
    if exitCode != 0:
      self.fail('combine_outputXX did not end correctly. Failing')

    ###########################
    ### CHECKING THE OUTPUT ###
    ###########################
    ds = nc.Dataset('outputs/schout_10.nc')
    # checking the time
    tmnc = ds.variables['time']
    lasttime = nc.num2date(tmnc[-1], tmnc.units, 'standard')
    self.assertEqual(datetime(2000, 1, 11), lasttime)

    xs = ds.variables['SCHISM_hgrid_node_x'][:]
    ys = ds.variables['SCHISM_hgrid_node_y'][:]

    ycond = np.logical_and(ys >= 33, ys <= 35)
    xcondFull = xs <= 6.5
    xcondAttenuated = xs >= 9
    
    # checking that hs is attenuated after the obstacles
    hs = ds.variables['WWM_1'][-1,:]
    self.assertTrue(np.all(hs[np.logical_and(ycond, xcondFull)] > 3.6))
    self.assertTrue(np.all(hs[np.logical_and(ycond, xcondAttenuated)] < 1))
   
    # checking that tm01 should not change
    tm = ds.variables['WWM_2'][-1,:]
    self.assertTrue(np.all(tm < 6.95))
    self.assertTrue(np.all(tm > 6.90))

    # closing the file
    ds.close()
    


    


if __name__ == '__main__':
  unittest.main()
