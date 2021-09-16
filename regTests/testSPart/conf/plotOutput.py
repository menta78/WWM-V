from netCDF4 import Dataset
from matplotlib import pyplot as plt
import numpy as np

from alphaBetaLab import abTriangularMesh as trmsh


def plotMapAtTime(ncfilepath, timeIndex):
  ds = Dataset(ncfilepath)
  hs = ds.variables['WWM_1'][timeIndex, :]
 #xs = ds.variables['SCHISM_hgrid_node_x'][:]
 #ys = ds.variables['SCHISM_hgrid_node_y'][:]
  msh = trmsh.loadFromGr3File('hgrid.gr3')
  ndid = list(msh.nodes.keys())
  ndid.sort()
  xs = np.array([msh.nodes[k][0] for k in ndid])
  ys = np.array([msh.nodes[k][1] for k in ndid])
  hs[hs.mask] = 0
  
  levels = np.arange(0., np.max(hs.flatten())*1.05, .02)
  cf = plt.tricontourf(xs, ys, hs, levels)
 #cf = bm.contourf(xs, ys, hs, tri=True)
  cf.cmap.set_over([.5,0,0])
  plt.triplot(xs, ys)
 #cf.set_clim(0, 1.6)
  plt.colorbar()
 #plt.savefig('plt_t=' + str(timeIndex) + '.png')
  plt.show()


if __name__ == '__main__':
  ncfilepath = 'outputs/schout_1.nc'
  timeIndex = -1
  import pdb; pdb.set_trace()
  plotMapAtTime(ncfilepath, timeIndex)

