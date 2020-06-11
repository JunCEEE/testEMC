import numpy as np
def logFactorial(n):
    if n < 20:
        value = np.log(np.math.factorial(n))
    else:
        value = 0.5*np.log(2*np.pi*n) + n*np.log(n/np.e)
    return value
def pnm(p,N,M):
    #this function gives the probability to observe a voxel at least M times
    #from an ensemble of N snapshots
    # p = probability to hit a voxel for a single snaprshot
    # N = Number of Snapshots
    # M = Redundancy
    if N < M:
        s = 1
    else:
        s = 0
        lp = np.log(p)
        lmp = np.log(1-p)
        for k in np.arange(M):
            s = s + np.exp(logFactorial(N) - logFactorial(N-k) - logFactorial(k) + k*lp + (N-k)*lmp)
    return np.maximum(1-s,0)
def numberOfSnapShots(d,D,nPhotons,SNR,P_tilde):
    #nS number of Snapshots
    #d is the Resolution
    # D is the Diameter of a Molecule
    #nPhotons is the Number of Photons
    #P_tilde is the combined probability
    #Number of Resolution elements
    R = D/d
    #number of voxels at Resolution Shell
    nV_Shell = 16*np.pi*R**2
    #probability per Shannon voxel
    p = 1./(4*R)
    M = np.ceil(SNR**2/nPhotons)
    # P -> Probability to observe a voxel at least M times from an
    # ensemble of nS snapshot
    # obtained from given P_tilde
    P = np.exp(2*np.log(P_tilde)/nV_Shell)

    print ("P =",P)
    print ("P_tilde =",P**(nV_Shell/2.0))
    print ("M =",M)
    print ("R =",R)
    print ("nV_Shell =",nV_Shell)

    #nSmax = 1e12#
    nSmax = 1e9#
    #step = 2**10#
    step = 100#
    nS0 = M
    while step > 1:
        for nS in np.arange(nS0,nSmax,step):
            if pnm(p,nS,M) > P:
                break
        nS0 = nS - step
        step = step /2
    return nS

if __name__ == "__main__":
    #ns = numberOfSnapShots(1,10,0.42e22,1,0.5)
    #ns = numberOfSnapShots(1,24.7,0.42e22,1,0.5)
    #ns = numberOfSnapShots(30,164,1.8e-2,1,0.5)
    #ns = numberOfSnapShots(25,164,6.7e-2,1,0.5)
    #ns = numberOfSnapShots(10,164,2.0e-3,1,0.5)
    #ns = numberOfSnapShots(25,164,6.7e-2,1,0.5)
    #ns = numberOfSnapShots(25,164,6.7e-2,1,0.5)
    #ns = numberOfSnapShots(30,164,1e-1,1,0.5)
    #ns = numberOfSnapShots(25,164,4.8e-2,1,0.5)
    #ns = numberOfSnapShots(10,164,2.6e-2,1,0.5)
    #ns = numberOfSnapShots(30,164,8.4e-3,1,0.5)
    #ns = numberOfSnapShots(25,164,3.8e-3,1,0.5)
    #ns = numberOfSnapShots(10,164,2.2e-3,1,0.5)
    #ns = numberOfSnapShots(25,164,4.8e-2,1,0.8)
    ns = numberOfSnapShots(10,164,2.6e-2,1,0.8)
    print (ns)
