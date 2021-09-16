import numpy as np
from datetime import datetime

from matplotlib import pyplot as plt

import testUtil


outNc = 'outputs/schout_1.nc'

pt = [6., 34.]
spfname = 'P-1.sp2d'
dateToLoad = datetime(2000, 1, 2)
timeIndex = -1

outpng = spfname + '.png'


p1hsPt, p1tmPt, p1dirmPt, p1dsprPt, p2hsPt, p2tmPt, p2dirmPt, p2dsprPt, p3hsPt, p3tmPt, p3dirmPt, p3dsprPt = testUtil.getPeaksAtPoint(outNc, timeIndex, pt)

tts, drs, spec = testUtil.getSpecAtDate(spfname, dateToLoad)

# correcting the direction for pcolor
#drs = drs + np.abs(np.diff(drs)[0])/2

def correctDir(pdir, drs):
  if pdir > np.max(drs):
    pdir = pdir-360
  if pdir < np.min(drs)-np.abs(np.diff(drs)[0]):
    pdir = pdir+360
  return pdir

p1dirmPt = correctDir(p1dirmPt, drs)
p2dirmPt = correctDir(p2dirmPt, drs)
p3dirmPt = correctDir(p3dirmPt, drs)

drs = np.concatenate((drs, np.array([drs[0]-360])), axis=0)
spec = np.concatenate((spec, np.reshape(spec[0,:], [1,25])), axis=0)

#fg = plt.figure()
#ax = fg.add_subplot(projection='polar')
#plt.contourf(drs/180*np.pi, tts, spec.transpose()**.5, 200, cmap='turbo')
#plt.ylim([1, 20])
#plt.scatter(np.array([p1dirmPt, p2dirmPt, p3dirmPt])/180*np.pi, [p1tmPt, p2tmPt, p3tmPt], s=[20*p1hsPt, 20*p2hsPt, 20*p3hsPt], c='red', edgecolors='lightgreen')

plt.style.use('dark_background')
fg = plt.figure()
ax = fg.add_subplot(projection='polar')
import pdb; pdb.set_trace()
fq = tts**-1
plt.contourf(drs/180*np.pi, fq, spec.transpose()**.5, 200, cmap='turbo')

ptfq = np.array([p1tmPt**-1, p2tmPt**-1 if p2tmPt != 0 else np.nan, p3tmPt**-1 if p3tmPt != 0 else np.nan])
maxfq = np.nanmax(ptfq)*1.3
plt.ylim([0, maxfq])

sprd = np.linspace(p1dirmPt-2*p1dsprPt, p1dirmPt+2*p1dsprPt, 300)/180*np.pi
pfrq = np.ones(sprd.shape)*p1tmPt**-1
plt.plot(sprd, pfrq, color='red', linewidth=1, zorder=1)

if p2tmPt != 0:
  sprd = np.linspace(p2dirmPt-2*p2dsprPt, p2dirmPt+2*p2dsprPt, 300)/180*np.pi
  pfrq = np.ones(sprd.shape)*p2tmPt**-1
  plt.plot(sprd, pfrq, color='red', linewidth=1, zorder=1)

if p3tmPt != 0:
  sprd = np.linspace(p3dirmPt-2*p3dsprPt, p3dirmPt+2*p3dsprPt, 300)/180*np.pi
  pfrq = np.ones(sprd.shape)*p3tmPt**-1
  plt.plot(sprd, pfrq, color='red', linewidth=1, zorder=1)

sc = plt.scatter(np.array([p1dirmPt, p2dirmPt, p3dirmPt])/180*np.pi, ptfq, 
        s=[20*p1hsPt, 20*p2hsPt, 20*p3hsPt], c='red', edgecolors='lightgreen', zorder=2)

plt.savefig(outpng, dpi=300)

plt.show()


