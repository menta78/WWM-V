import numpy as np
from scipy.io import FortranFile
from datetime import datetime

import xarray as xr
from matplotlib.tri import Triangulation, LinearTriInterpolator



def getSpecAtDate(specBinFileName, dateToLoad):
  dateToLoadStr = dateToLoad.strftime('%Y%m%d.%H%M%S')

  rkind = 'float64'
  
  fl = FortranFile(specBinFileName, 'r')
  numsig = fl.read_ints()[0]
  numdir = fl.read_ints()[0]
  
  spsig = fl.read_reals(dtype=rkind)
  tm = 2*np.pi/spsig
  spdir = fl.read_reals(dtype=rkind)
  
  isum = fl.read_ints()
  
  datestr = ''
  while datestr != dateToLoadStr:
    dateAsByte = fl.read_record(dtype='byte')
    datestr = dateAsByte.tostring().decode('ascii')
    
    dep = fl.read_reals(dtype=rkind)
    curtxy = fl.read_reals(dtype=rkind)
    spec = np.reshape(fl.read_reals(dtype=rkind), [numdir, numsig])
    acloc = np.reshape(fl.read_reals(dtype=rkind), [numdir, numsig])

  fl.close()

  frq = spsig/2/np.pi 
  tts = frq**-1
  drs = (np.pi/2-spdir)/np.pi*180+180

  return tts, drs, spec



def getPeaksAtPoint(wwmOutNc, timeIndex, pt):
  ds = xr.open_dataset(wwmOutNc)
  p1hs = ds['WWM_P1HS'][timeIndex,:].values
  p2hs = ds['WWM_P2HS'][timeIndex,:].values
  p3hs = ds['WWM_P3HS'][timeIndex,:].values
  p1tm = ds['WWM_P1TM01'][timeIndex,:].values
  p2tm = ds['WWM_P2TM01'][timeIndex,:].values
  p3tm = ds['WWM_P3TM01'][timeIndex,:].values
  p1dirm = ds['WWM_P1DM'][timeIndex,:].values
  p2dirm = ds['WWM_P2DM'][timeIndex,:].values
  p3dirm = ds['WWM_P3DM'][timeIndex,:].values
  p1dspr = ds['WWM_P1DSPR'][timeIndex,:].values
  p2dspr = ds['WWM_P2DSPR'][timeIndex,:].values
  p3dspr = ds['WWM_P3DSPR'][timeIndex,:].values
  xs = ds['SCHISM_hgrid_node_x'][:].values
  ys = ds['SCHISM_hgrid_node_y'][:].values
  ds.close()

  triObj = Triangulation(xs,ys)
  intpltr = LinearTriInterpolator(triObj, p1hs)
  p1hsPt = intpltr(pt[0], pt[1])
  intpltr = LinearTriInterpolator(triObj, p1tm)
  p1tmPt = intpltr(pt[0], pt[1])
  intpltr = LinearTriInterpolator(triObj, p1dirm)
  p1dirmPt = intpltr(pt[0], pt[1])
  intpltr = LinearTriInterpolator(triObj, p1dspr)
  p1dsprPt = intpltr(pt[0], pt[1])

  intpltr = LinearTriInterpolator(triObj, p2hs)
  p2hsPt = intpltr(pt[0], pt[1])
  intpltr = LinearTriInterpolator(triObj, p2tm)
  p2tmPt = intpltr(pt[0], pt[1])
  intpltr = LinearTriInterpolator(triObj, p2dirm)
  p2dirmPt = intpltr(pt[0], pt[1])
  intpltr = LinearTriInterpolator(triObj, p2dspr)
  p2dsprPt = intpltr(pt[0], pt[1])

  intpltr = LinearTriInterpolator(triObj, p3hs)
  p3hsPt = intpltr(pt[0], pt[1])
  intpltr = LinearTriInterpolator(triObj, p3tm)
  p3tmPt = intpltr(pt[0], pt[1])
  intpltr = LinearTriInterpolator(triObj, p3dirm)
  p3dirmPt = intpltr(pt[0], pt[1])
  intpltr = LinearTriInterpolator(triObj, p3dspr)
  p3dsprPt = intpltr(pt[0], pt[1])

  return float(p1hsPt), float(p1tmPt), float(p1dirmPt), float(p1dsprPt), float(p2hsPt), float(p2tmPt), float(p2dirmPt), float(p2dsprPt), float(p3hsPt), float(p3tmPt), float(p3dirmPt), float(p3dsprPt)
