import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import svds

edges = np.loadtxt('edges.txt', dtype=int)
rows = edges[:, 1] - 1 
cols = edges[:, 0] - 1
data = np.ones(len(edges))

A = coo_matrix((data, (rows, cols)), shape=(8297, 8297))

u, s, vt = svds(A.asfptype(), k=1)

auths = u.flatten()
hubs = vt.flatten()

if np.sum(hubs) < 0: hubs *= -1
if np.sum(auths) < 0: auths *= -1

with open('svd_results.txt', 'w') as f:
    f.write(' '.join(map(str, auths)) + '\n')
    f.write(' '.join(map(str, hubs)) + '\n')
    f.write(str(s[0]) + '\n')