import os
import numpy as np
from matplotlib import pyplot as plt

# importing from alphaBetaLab the needed components
from alphaBetaLab import abTriangularMesh, abHighResAlphaMatrix, abTriangularMeshGridBuilder, abEstimateAndSave
from alphaBetaLab.abOptionManager import abOptions


# definition of the spectral grid
dirs = np.linspace(0, 2*np.pi, 9)
nfreq = 5
minfrq = .1
frqfactor = 1.1
freqs = [minfrq*(frqfactor**i) for i in range(1,nfreq + 1)]

# definition of the spatial mesh
gridname = 'test'
mshfile = 'hgrid.gr3'
triMeshSpec = abTriangularMesh.loadFromGr3File(mshfile)

# output directory
outputDestDir = './'

# number of cores for parallel computing
nParWorker = 12
nParWorker = 2

# this option indicates that the computation should be skipped for cells smaller than 3 km
minSizeKm = 3
opt = abOptions(minSizeKm=minSizeKm)

# create fake bathymetry for the test
meshxs = [xy[0] for xy in triMeshSpec.nodes.values()]
meshys = [xy[1] for xy in triMeshSpec.nodes.values()]
x = np.linspace(np.min(meshxs), np.max(meshxs), 50)
y = np.linspace(np.min(meshys), np.max(meshys), 50)
z = -np.ones((len(y), len(x)))*500
# adding the obstructions
obstIndxs = np.arange(25)*2
z[obstIndxs, 15] = 0
alphamtx = np.ones(z.shape)
alphamtx[z > -.1] = 0
plt.pcolor(x, y, alphamtx)
plt.triplot(meshxs, meshys)
plt.show()

highResolutionBathyMatrix = abHighResAlphaMatrix.abHighResAlphaMatrix(x, y, alphamtx)

# create mesh builder
gridBld = abTriangularMeshGridBuilder.abTriangularMeshGridBuilder(triMeshSpec, nParWorker = nParWorker)
grid = gridBld.buildGrid()


# instruction to do the computation and save the output
abEstimateAndSave._abEstimateAndSave(dirs, freqs, gridname, grid, highResolutionBathyMatrix, outputDestDir, nParWorker, abOptions=opt)
os.system('mv obstructions_local.test.in obstructions_local.in; mv obstructions_shadow.test.in obstructions_shadow.in')




