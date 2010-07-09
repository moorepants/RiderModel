from RiderLeanAB import RiderLeanAB
import matplotlib.pyplot as plt
import numpy as np

v = np.linspace(0., 10., num=100)

for speed in v:
    A, B, bike = RiderLeanAB(8, speed)
    w = np.linalg.eig(A)[0]
    #print "Computed the eigenvalues at", speed
    #print "Eigenvalues\n", w
    try:
        eigvals = np.vstack((eigvals, w))
    except:
        eigvals = w

plt.plot(v, eigvals, 'k.')
plt.ylim(-10,10)
plt.title(bike)
plt.show()
