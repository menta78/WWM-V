# author: Lorenzo Mentaschi
import os
import unittest
import netCDF4 as nc
import numpy as np
from datetime import datetime

from wwmTestUtil import wwmTestUtil as util

from testSPart import testUtil


class testSPart(util.wwmTestTemplate):

  def test(self):
    # launching schismWWM
    launchCmd = util.getScWwmLaunchCommand(nProcesses=2)
    exitCode = os.system(launchCmd)
    if exitCode != 0:
      self.fail('schismWWM did not end correctly. Failing')
    # combining the output
    combineCmd = util.getCombineCommand(bgnParam=1, endParam=1)
    exitCode = os.system(combineCmd)
    if exitCode != 0:
      self.fail('combine_outputXX did not end correctly. Failing')

    ###########################
    ### CHECKING THE OUTPUT ###
    ###########################
    ds = nc.Dataset('outputs/schout_1.nc')

    hs = ds.variables['WWM_1'][:]

    # checking that all the fields have been saved
    p1hs = ds.variables['WWM_P1HS'][:]
    p1tm01 = ds.variables['WWM_P1TM01'][:]
    p1dm = ds.variables['WWM_P1DM'][:]
    p1dspr = ds.variables['WWM_P1DSPR'][:]

    p2hs = ds.variables['WWM_P2HS'][:]
    p2tm01 = ds.variables['WWM_P2TM01'][:]
    p2dm = ds.variables['WWM_P2DM'][:]
    p2dspr = ds.variables['WWM_P2DSPR'][:]

    p3hs = ds.variables['WWM_P3HS'][:]
    p3tm01 = ds.variables['WWM_P3TM01'][:]
    p3dm = ds.variables['WWM_P3DM'][:]
    p3dspr = ds.variables['WWM_P3DSPR'][:]

    # closing the file
    ds.close()

    self.assertTrue(np.all(p1hs**2 + p2hs**2 + p3hs**2 - hs**2 <= 1e-5), 'the peaks\'s energy does not sum up to hs')
    self.assertTrue(np.all(p1hs >= p2hs), 'some problem with peak ordering')
    self.assertTrue(np.all(p2hs >= p3hs), 'some problem with peak ordering')
    ratioP1NotZero = np.sum(p1hs>0)/len(p1hs.flatten())
    self.assertTrue(ratioP1NotZero > .999)
    ratioP2NotZero = np.sum(p2hs>0)/len(p2hs.flatten())
    self.assertTrue(ratioP2NotZero > .8) # the second peak develops after some fetch
    ratioP3NotZero = np.sum(p3hs>0)/len(p3hs.flatten())
    self.assertTrue(ratioP3NotZero < .03) # p3hs should be always 0 but in about 2% of the point/time the algorithm splits one peak

    pt = [6., 34.]
    timeIndex = -1
    # checking the value at the location where the spectrum was saved
    p1hsPt, p1tmPt, p1dirmPt, p1dsprPt, p2hsPt, p2tmPt, p2dirmPt,\
    p2dsprPt, p3hsPt, p3tmPt, p3dirmPt, p3dsprPt =\
                     testUtil.getPeaksAtPoint('outputs/schout_1.nc', timeIndex, pt)

    self.assertTrue(p1hsPt > 1.8)
    self.assertTrue(p1tmPt > 4.5)
    self.assertTrue(160<p1dirmPt and p1dirmPt<180)
    self.assertTrue(p1dsprPt > 40)

    self.assertTrue(p2hsPt < 1)
    self.assertTrue(p2tmPt > 11)
    self.assertTrue(260<p2dirmPt and p1dirmPt<280)
    self.assertTrue(p2dsprPt < 30)

    self.assertTrue(p3hsPt == 0)
    self.assertTrue(p3tmPt == 0)
    self.assertTrue(p3dirmPt == 0)
    self.assertTrue(p3dsprPt == 0)

    
    

    


    


if __name__ == '__main__':
  unittest.main()
