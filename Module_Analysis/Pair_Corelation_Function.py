import numpy as np

def compute_2_body_dist_prob(R,bins=100):
    maxR = np.max(R)
    ddist = maxR/(bins-1)
    def I(dist):
        return int(dist/ddist)
    def v_shell(dist):
        return 4/3*np.pi*((dist+ddist)**3 - dist**3)
    X = np.linspace(0,maxR,bins)
    PR = np.zeros(bins,dtype=float)
    for r in R:
        for r1,r2 in zip(r[:r.shape[0]-1],r[1:]):
            if any(r1!=r2):
                try:
                    PR[I(np.linalg.norm(r1-r2))] +=1/v_shell(np.linalg.norm(r1-r2))*1/(R.shape[0]*R.shape[1])
                except IndexError:
                    pass
                    #PR.resize(I(np.linalg.norm(r1-r2))+1)
                    #PR[I(np.linalg.norm(r1-r2))] +=1/v_shell(np.linalg.norm(r1-r2))*1/(R.shape[0]*R.shape[1])
                    #X = np.linspace(0,maxR,I(np.linalg.norm(r1-r2))+1)
    return X,PR
def compute_2_body_dist_prob_from_dist(D,bins=100):
    maxD = np.max(D)
    ddist = maxD/(bins-1)
    def I(dist):
        return int(dist/ddist)
    def v_shell(dist):
        return 4/3*np.pi*((dist+ddist)**3 - dist**3)
    X = np.linspace(0,maxD,bins)
    PR = np.zeros(bins,dtype=float)
    for D in D:
        if D!=0:
            try:
                PR[I(D)] += 1/v_shell(D)*1/(D.shape[0]*D.shape[1])
            except IndexError:
                pass
    return X,PR