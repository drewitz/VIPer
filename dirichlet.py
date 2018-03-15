#!/usr/bin/env python3

# should draw the dirichlet fundamental domain

import uppergeodesic as geod
import numpy as np
import numpy.linalg as lina
import matplotlib.pyplot as plt
import mpmath as mp

exp = mp.exp
log = mp.log
sqrt = mp.sqrt


def dirichlet(p,matrices):
    """
    takes p as the point and applies matrices.
    Then draws p, it's images under the matrices and the bisectors
    """
    geodesics = []
    for A in matrices:
        geodesics.append(get_bisector(p, mob(A,p)))


def get_bisector(a,b):
    """
    a and b are tuples (real, imag)
    """
    assert a!=b, "no bisector for a point and itself"

    if a[0] == b[0]:
        return geod.UpperGeodesic.from_midpoint_and_radius(a[0],
                a[1]*exp(1/2*log(b[1]/a[1]))
                )
    elif a[1] == b[1]:
        return geod.UpperGeodesic.vertical((b[0]+a[0])/2)
    else:
        impart = a[1]*exp(1/2*log(b[1]/a[1]))
        realpart = (b[0]-a[0])/(b[1]-a[1])*b[1]
        m = b[0] - realpart
        r = sqrt(impart**2 + realpart**2)
        return geod.UpperGeodesic.from_midpoint_and_radius(m,r)

def mob(A,p):
    """
    matrix A acts on p
    """
    pcomp = p[0]+ p[1]*1j
    tcomp = (A[0,0]*pcomp + A[0,1]) / (A[1,0]*pcomp + A[1,1])
    return (tcomp.real, tcomp.imag)



#firstmatrices = [np.matrix([[1,2],[0,1]]), np.matrix([[1,0],[2,1]]),
#        np.matrix([[1,2],[2,5]]), np.matrix([[1,-2],[-2,5]]),
#        np.matrix([[1,-2],[2,-3]]), np.matrix([[1,2],[-2,-3]]),
#        np.matrix([[5,2],[2,1]]), np.matrix([[5,-2],[-2,1]]),
#        np.matrix([[-3,-2],[2,1]]), np.matrix([[-3,2],[-2,1]]),
#        np.matrix([[3,2],[4,3]]), np.matrix([[3,4],[2,3]])]

#matrices = firstmatrices + [lina.inv(A) for A in firstmatrices]

matrices = [np.matrix([[1,2],[0,1]]), np.matrix([[1,0],[2,1]]),
            np.matrix([[1,-2],[0,1]]), np.matrix([[1,0],[-2,1]])]

dirichlet((0.5,2), matrices)

geod.UpperGeodesic.plot_all()
plt.show()
