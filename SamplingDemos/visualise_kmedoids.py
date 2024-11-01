import numpy as np
import matplotlib.pyplot as plt

from sklearn_extra.cluster import KMedoids
from quippy.clustering import get_cur_scores
from skmatter.feature_selection import CUR, FPS
from scipy.spatial.distance import cdist
from scipy.spatial import ConvexHull

npoints = 31

points = np.zeros((npoints * npoints, 2))

idx = 0
for i in range(npoints):
    for j in range(npoints):
        points[idx, :] = [i/(npoints-1), j/(npoints-1)]
        idx += 1




N = 20

colors = [f"C{i}" for i in range(10)]
markers = [".", "o", "*"]

colmar = [(c, m) for c in colors for m in markers]

np.random.seed(42)

method = "KMED"

for kind in ["uniform", "quadratic"]:
    plt.clf()
    
    km = KMedoids(n_clusters=N)
    km.fit(points)
    x = points[km.medoid_indices_, :]

    idxs = km.predict(points)
    
    for i in range(N):
        col, mar = colmar[i]
        plt.scatter(points[idxs==i, 0], points[idxs==i, 1], color=col, marker=mar)

            
    dists = cdist(x, x).flatten()
    av_sep = np.average(dists[dists>0])

    hull = ConvexHull(x, False)

    plt.scatter(x[:, 0], x[:, 1], marker="x", color=f"k")
    plt.title(method + f"; Coverage: {100*hull.volume:.1f}%, Average Separation: {np.round(av_sep, 3)}")
    plt.savefig(f"{method}_{kind}.png", dpi=400)
    
    points = points ** 2