from alphaBetaLab import abTriangularMesh as trmsh
from matplotlib import pyplot as plt

msh = trmsh.loadFromGr3File('hgrid.gr3')
xs = [msh.nodes[k][0] for k in msh.nodes.keys()]
ys = [msh.nodes[k][1] for k in msh.nodes.keys()]

plt.triplot(xs, ys, linewidth=.1, color='k')


bcmsh = trmsh.loadFromGr3File('wwmbnd.gr3')
xsbdc = [bcmsh.nodes[k][0] for k in bcmsh.nodes.keys()]
ysbdc = [bcmsh.nodes[k][1] for k in bcmsh.nodes.keys()]
bt = [bcmsh.nodeBathy[k] for k in bcmsh.nodes.keys()]

plt.scatter(xsbdc, ysbdc, c=bt)

plt.show()

