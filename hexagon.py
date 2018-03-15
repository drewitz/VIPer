#!/usr/bin/env python3

import uppergeodesic as gd
import numpy as np

# should plot a hexagon in upper half space given three parameters whose
# logarithms are the side lengths
d1 = -0.3
d2 = 3
d3 = 0.3

# UpperGeodesic doesn't keep track of the orientation!!! :(
s1 = gd.UpperGeodesic(-1,1, label='S1, q1={}'.format(d1))
s2 = gd.UpperGeodesic(0, 'inf', 'r', label='S2, q2={}'.format(d2))
s3 = gd.UpperGeodesic(-d2, d2, 'g', label='S3, q3={}'.format(d3))

temp4 = (d3 - 1)/(d3 + 1)
s4a = d2*temp4
s4b = d2/temp4
s4 = gd.UpperGeodesic(s4a, s4b, 'k', label = 'S4')

s6a = (1 - d1)/(1 + d1)
s6b = 1/s6a

# 
s5a = (s4a*s4b*s6a + s6a*s6b**2 - (np.sqrt(-s4a + s6a)*np.sqrt(-s4b + s6a)*s6a + np.sqrt(-s4a + s6a)*np.sqrt(-s4b + s6a)*s6b)*np.sqrt(-s4a + s6b)*np.sqrt(-s4b + s6b) + (s4a*s4b - 2*(s4a + s4b)*s6a + s6a**2)*s6b)/(2*s4a*s4b - 2*np.sqrt(-s4a + s6a)*np.sqrt(-s4a + s6b)*np.sqrt(-s4b + s6a)*np.sqrt(-s4b + s6b) - (s4a + s4b)*s6a + s6a**2 - (s4a + s4b)*s6b + s6b**2)
s5b = (s4a*s4b - s6a*s6b - np.sqrt(s4a**2*s4b**2 + s4a*s4b*s6a**2 + (s4a*s4b - (s4a + s4b)*s6a + s6a**2)*s6b**2 - (s4a**2*s4b + s4a*s4b**2)*s6a - (s4a**2*s4b + s4a*s4b**2 + (s4a + s4b)*s6a**2 - (s4a**2 + 2*s4a*s4b + s4b**2)*s6a)*s6b))/(s4a + s4b - s6a - s6b)

s5 = gd.UpperGeodesic(s5a, s5b, 'm', label = 'S5')
s6 = gd.UpperGeodesic(s6a, s6b, 'y', label = 'S6')

gd.UpperGeodesic.plot_all()
gd.plt.show()
