import numpy as np
import netCDF4 as nc

from alphaBetaLab import abTriangularMesh as tri

outFlName = 'ice.nc'
lonlat2d = False

msh = tri.loadFromGr3File('hgrid.gr3')
xs = np.array([xy[1][0] for xy in msh.nodes.items()])
ys = np.array([xy[1][1] for xy in msh.nodes.items()])

minx, maxx = np.min(xs), np.max(xs)
miny, maxy = np.min(ys), np.max(ys)

xsice = np.arange(minx-.1, maxx+.2, .1)
ysice = np.arange(miny-.1, maxy+.2, .1)

xsmtx, ysmtx = np.meshgrid(xsice, ysice)

fl = nc.Dataset(outFlName, 'w')
lat_dim = fl.createDimension('y', len(ysice)) # latitude axis
lon_dim = fl.createDimension('x', len(xsice)) # longitude axis
time_dim = fl.createDimension('time', None) # unlimited axis (can be appended to).
fl.title = 'test ice data'

if lonlat2d:
  lat = fl.createVariable('lat', np.float32, ('y', 'x'))
  lat.units = 'degrees_north'
  lat.long_name = 'latitude'
  lat[:] = ysmtx
  
  lon = fl.createVariable('lon', np.float32, ('y', 'x'))
  lon.units = 'degrees_north'
  lon.long_name = 'longitude'
  lon[:] = xsmtx
else:
  lat = fl.createVariable('lat', np.float32, ('y'))
  lat.units = 'degrees_north'
  lat.long_name = 'latitude'
  lat[:] = ysice
  
  lon = fl.createVariable('lon', np.float32, ('x'))
  lon.units = 'degrees_north'
  lon.long_name = 'longitude'
  lon[:] = xsice

tmnc = fl.createVariable('time', np.float32, ('time',))
tmnc.units = 'days since 2000-01-01'
tmnc.long_name = 'time'
tms = np.arange(0, 6, 1) #5 days of ice every 24 hours
tmnc[:] = tms

icenc = fl.createVariable('siconc', np.float32, ('time', 'y', 'x'))
icenc.units = 'surface ratio'
icenc.long_name = 'sea ice concentration'
icenc.coordinates = 'lon lat'

offsetDeg = .5
factor = .5
for itm in range(len(tms)):
  icec0 = (xsice-minx-offsetDeg*itm)*factor
  icec0[icec0<0] = 0
  icec0[icec0>1] = 1
  iceconc = np.tile(icec0, [len(ysice), 1])
  icenc[itm, :] = iceconc

fl.close()


