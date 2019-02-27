import numpy as np


def edge_topology(F):
    uids = {}
    EV = []
    EF = []
    d = F.shape[1]
    FE = np.ones(shape=(F.shape[0], d), dtype=np.int32) * -1
    nEdges = 0

    for i in range(0, F.shape[0]):
        f = F[i]
        vtx = f[f >= 0]
        n = len(vtx)
        for j in range(0, n):
            e = (vtx[j], vtx[(j+1) % n])
            if e[0] < e[1]:
                uid = (e[0], e[1])
                side = 0
            else:
                uid = (e[1], e[0])
                side = 1
            if uid not in uids:
                uids[uid] = nEdges
                EV.append(e)
                ef = [-1, -1]
                ef[side] = i
                EF.append(ef)
                FE[i, j] = nEdges
                nEdges += 1
            else:
                eid = uids[uid]
                if EF[eid][side] != -1:
                    raise Exception('Bad Triangle Oriantation triangles %i and %i' %(i, EF[eid][side]))
                EF[eid][side] = i
                FE[i, j] = eid

    return np.array(EV, dtype=np.int32), np.array(FE, dtype=np.int32), np.array(EF, dtype=np.int32)

