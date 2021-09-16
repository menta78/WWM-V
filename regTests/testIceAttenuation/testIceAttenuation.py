# author: Lorenzo Mentaschi
import os
import unittest
import netCDF4 as nc
import numpy as np
from datetime import datetime

from wwmTestUtil import wwmTestUtil as util


class testIceAttenuation(util.wwmTestTemplate):

  def test(self):
    # launching schismWWM
    launchCmd = util.getScWwmLaunchCommand(nProcesses=2)
    exitCode = os.system(launchCmd)
    if exitCode != 0:
      self.fail('schismWWM did not end correctly. Failing')
    # combining the output for the first 3 hours of the first day
    combineCmd = util.getCombineCommand(bgnParam=1, endParam=1)
    exitCode = os.system(combineCmd)
    if exitCode != 0:
      self.fail('combine_outputXX did not end correctly. Failing')
    # combining the output for the last 3 hours of the last day
    combineCmd = util.getCombineCommand(bgnParam=8, endParam=8)
    exitCode = os.system(combineCmd)
    if exitCode != 0:
      self.fail('combine_outputXX did not end correctly. Failing')

    ###########################
    ### CHECKING THE OUTPUT ###
    ###########################

    # 1st file
    ds = nc.Dataset('outputs/schout_1.nc')
    # checking the time
    tmnc = ds.variables['time']
    lasttime = nc.num2date(tmnc[2], tmnc.units, 'standard')
    self.assertEqual(datetime(2000, 1, 1, 3, 0, 0), lasttime)
    
    # checking that ice is loaded correctly and hs is duely attenuated
    icec = ds.variables['WWM_ice_conc'][-1,:]
    hs = ds.variables['WWM_1'][-1,:]
    xs = ds.variables['SCHISM_hgrid_node_x'][:]
    ys = ds.variables['SCHISM_hgrid_node_y'][:]
    ycond = np.logical_and(ys >= 33, ys <= 35)
    xcondNoIce = xs <= 4.5
    xcondAttenuated = xs >= 7
    # icec
    self.assertTrue(np.all(icec[np.logical_and(ycond, xcondNoIce)] <= .000001))
    self.assertTrue(np.all(icec[np.logical_and(ycond, xcondAttenuated)] > .999))
    # hs
    self.assertTrue(np.all(hs[np.logical_and(ycond, xcondNoIce)] >= 3.5))
    self.assertTrue(np.all(hs[np.logical_and(ycond, xcondAttenuated)] < .3))
    ds.close()
    
    # last file
    ds = nc.Dataset('outputs/schout_8.nc')
    # checking the time
    tmnc = ds.variables['time']
    lasttime = nc.num2date(tmnc[-1], tmnc.units, 'standard')
    self.assertEqual(datetime(2000, 1, 5, 0, 0, 0), lasttime)
    
    # checking that ice is loaded correctly and hs is duely attenuated
    icec = ds.variables['WWM_ice_conc'][-1,:]
    hs = ds.variables['WWM_1'][-1,:]
    xs = ds.variables['SCHISM_hgrid_node_x'][:]
    ys = ds.variables['SCHISM_hgrid_node_y'][:]
    ycond = np.logical_and(ys >= 33, ys <= 35)
    xcondNoIce = xs <= 5.5
    xcondAttenuated = xs >= 9.5
    # icec
    self.assertTrue(np.all(icec[np.logical_and(ycond, xcondNoIce)] <= .001))
    self.assertTrue(np.all(icec[np.logical_and(ycond, xcondAttenuated)] > .999))
    # hs
    self.assertTrue(np.all(hs[np.logical_and(ycond, xcondNoIce)] >= 3.5))
    self.assertTrue(np.all(hs[np.logical_and(ycond, xcondAttenuated)] < .1))
    ds.close()
    


    


if __name__ == '__main__':
  unittest.main()
