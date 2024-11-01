import numpy as np
import matplotlib.pyplot as plt

from sklearn_extra.cluster import KMedoids
from quippy.clustering import get_cur_scores
from skmatter.feature_selection import CUR, FPS
from scipy.spatial.distance import cdist
from scipy.spatial import ConvexHull

def kmed(x, N):
    kmed = KMedoids(n_clusters=N, method="pam", init="k-medoids++").fit(x)
    return x[kmed.medoid_indices_, :]

def fps(x, N):
    selector = FPS(n_to_select=N)
    return selector.fit_transform(x.T).T

def cur(x, N):
    scores = get_cur_scores(x, clip_scores=False)

    idxs = np.arange(x.shape[0])
    sidx = np.random.choice(idxs, N, replace=False, p=scores/np.sum(scores))
    return x[sidx, :]
    

npoints = 31

points = np.zeros((npoints * npoints, 2))

idx = 0
for i in range(npoints):
    for j in range(npoints):
        points[idx, :] = [i/(npoints-1), j/(npoints-1)]
        idx += 1

methods = {
    "FPS" : fps,
    "CUR" : cur
}


markers = ["x", "*", "^" , "s"]

N = 20

np.random.seed(42)

for kind in ["uniform", "quadratic"]:
    for method in methods.keys():
        print(method)
        plt.clf()
        plt.scatter(points[:, 0], points[:, 1], marker=".")

        n_runs = 1 if method in ["KMED", "FPS"] else 4

        av_seps = []
        V = []

        for i in range(n_runs):
            x = methods[method](points, N)
            
            dists = cdist(x, x).flatten()
            av_sep = np.average(dists[dists>0])
            av_seps.append(av_sep)

            hull = ConvexHull(x, False)
            V.append(hull.volume)

            plt.scatter(x[:, 0], x[:, 1], marker=markers[i])
        plt.title(method + f"; Coverage: {100*np.average(V):.1f}%, Average Separation: {np.round(np.average(av_seps), 3)}")
        plt.savefig(f"{method}_{kind}.png", dpi=400)
    
    points = points ** 2