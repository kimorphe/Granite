from scipy.spatial import Delaunay, delaunay_plot_2d, Voronoi, voronoi_plot_2d
import matplotlib.pyplot as plt
import numpy as np

w=h=360

n=20

np.random.seed(0)
pts=np.random.randint(0,w,(n,2))

print(pts)

tri=Delaunay(pts)
fig=delaunay_plot_2d(tri)

vor=Voronoi(pts)
fig=voronoi_plot_2d(vor)

print(vor.vertices)
print(vor.regions)
plt.show()
