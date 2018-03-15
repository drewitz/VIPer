#!/usr/bin/env python3

import uppergeodesic as gd
import numpy as np
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('lengths', nargs = 2, type=float)
args = parser.parse_args()
# try 0.3 and 3

# should plot a pentagon in upper half space given two parameters whose
# logarithms are the side lengths

d1 = args.lengths[0]
d2 = args.lengths[1]

print(d1)
print(d2)

# UpperGeodesic doesn't keep track of the orientation!!! :(
s1 = gd.UpperGeodesic('inf', 0, label='S1, q1 = {0}'.format(d1))
s2 = gd.UpperGeodesic(-1, 1, 'r', label='S2, q2 = {0}'.format(d2))
s3b = (1-d2) / (1+d2)
s3a = 1/s3b
print(s3b)
s3 = gd.UpperGeodesic(s3a, s3b, 'y', label='S3')
s5a = -d1
s5b = d1

p = (1 + s5a**2)/(s3a+s3b)
q = s5a**2

temp = p + np.sqrt(p**2 - q)

s4 = gd.UpperGeodesic(s5a**2/temp, temp, 'm', label='S4')
s5 = gd.UpperGeodesic(s5a, s5b, 'g', label='S5')

gd.UpperGeodesic.plot_all()
gd.plt.show()
