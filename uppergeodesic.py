#!/usr/bin/env python3

# class file uppergeodesic.py

# started as a script to visualize what happens to hyperbolic plane
# if different isometries act on it

import geodesic as gd
import numpy as np
import numpy.linalg as lina
import matplotlib.pyplot as plt

# upper half space as the basic model

class UpperGeodesic(gd.Geodesic):
    """UpperGeodesic line in upper half space

    takes endpoints on boundary as arguments
    stores x and y data as points in x and y
    string "inf" is point at infinity, i.e. y=inf
    """
    
    xmin = 0
    xmax = 0
    ymin = 0 # just for consistency, shouldn't change
    ymax = 0
    inf = "inf"

    def __init__(self, a, b, color="b", label=''):
        """initialize UpperGeodesic by endpoints

        a, b - x values of the endpoints or "inf" if infinity
        res is resolution
        """

        super().__init__(a, b, color, label)
        # adjust the boundaries of hyperbolic space
        if self.start != UpperGeodesic.inf:
            if self.start < UpperGeodesic.xmin:
                UpperGeodesic.xmin = self.start
        if self.end > UpperGeodesic.xmax:
            UpperGeodesic.xmax = self.end
        UpperGeodesic.ymax = (UpperGeodesic.xmax - UpperGeodesic.xmin)/2

    @classmethod
    def vertical(cls, real):
        return cls(cls.inf, real)

    @classmethod
    def from_midpoint_and_radius(cls, m, r):
        """
        m is only the real part of the circle thing
        """
        return cls(m-r, m+r)

    def sort_se(self):
        """sort start and end"""
        if self.end == self.inf:
            # just want to assume that the first value is inf if any
            self.end = self.start
            self.start = self.inf
        if self.start != self.inf and self.end < self.start:
            # swap a and self.end such that a < self.end
            c = self.start
            self.start = self.end
            self.end = c

    def get_data(self):
        if self.start == UpperGeodesic.inf:
            # vertical line
            xs = [self.end, self.end]
            ys = [self.ymin, self.ymax]
        else:
            # calculate semicircle
            t = np.linspace(0, np.pi, self.res)
            r = (self.end - self.start)/2
            xs = r*(1 + np.cos(t)) + self.start
            ys = r*np.sin(t)
        return(xs, ys)

    ## the next two functions create new geodesics from existing ones
    def new_geod(self, a, b, c, d):
        """return new geodesic by received by moebius trafo

        apply the matrix
        | a b |
        | c d |
        on the geodesic self and return the resulting geodesic
        """

        start = self.apply_moebius(a, b, c, d, self.start)
        end = self.apply_moebius(a, b, c, d, self.end)

        return(UpperGeodesic(start, end))

    def new_from_matrix(self, M):
        return self.new_geod(M[0,0], M[0,1], M[1,0], M[1,1])

    ## apply transformations to ONE geodesic
    def apply_matrix(self, M):
        self.start = self.apply_moebius(M[0,0], M[0,1], M[1, 0], M[1,1],
                                        self.start)
        self.end = self.apply_moebius(M[0,0], M[0,1], M[1, 0], M[1,1],
                                      self.end)
        self.sort_se()

    def translate_one_geod(self, dx):
        if self.start != UpperGeodesic.inf:
            self.start += dx
        if self.end != UpperGeodesic.inf:
            self.end += dx

    def translate_one_at_zero(self, dx):
        """inverts at unit sphere, translates and inverts again"""
        a = self.inversion_on_unit_circle(self.start)
        b = self.inversion_on_unit_circle(self.end)
        if a != UpperGeodesic.inf:
            a += dx
        if b != UpperGeodesic.inf:
            b += dx
        self.start = self.inversion_on_unit_circle(a)
        self.end = self.inversion_on_unit_circle(b)
        self.sort_se()

    def rotate_one_geod(self, phi):
        """rotates the geodesic on upper half space
        
        conjugate to a rotation around origin in the disc model
        """
        if self.start == UpperGeodesic.inf:
            alpha = -np.pi/2
        else:
            alpha = self.from_upper_to_disc(self.start)
        beta = self.from_upper_to_disc(self.end)

        alpha += phi
        beta += phi

        self.start = self.from_disc_to_upper(alpha)
        self.end = self.from_disc_to_upper(beta)
        self.sort_se()

    def hyperbolic_translate_one(self, dmult=1.001):
        """translates one geodesic along UpperGeodesic(-1,1)"""
        diag = (dmult + 1.0/dmult)/2.0
        off = (1.0/dmult - dmult)/2.0
        matrix = np.matrix([[diag, off], [off, diag]])
        self.apply_matrix(matrix)

    # tesselate hyperbolic space
    @classmethod
    def tesselate(self, depth=10):
        """Tesselates according to SL(2,Z)"""
        g0 = UpperGeodesic(-1,1, "r")
        g1 = UpperGeodesic(-0.5,self.inf, "r")
        g2 = UpperGeodesic(0.5,self.inf, "r")
        first = [g0,g1,g2]
        for k in range(1, depth):
            for g in first:
                g.new_geod(1, k, 0, 1)
                g.new_geod(1, -k, 0, 1)
        kmax = len(UpperGeodesic.all_geods)
        for geod in UpperGeodesic.all_geods[:kmax]:
            temp = [geod.new_geod(0, -1, 1, 0)]
            for k in range(1, 2*depth):
                temp.append(geod.new_geod(1, 0, k, 1))
                temp.append(geod.new_geod(1, 0, -k, 1))
            for k in range(1, depth//2):
                for t in temp:
                    t.new_geod(1, k, 0, 1)
                    t.new_geod(1, -k, 0, 1)

        UpperGeodesic.xmin= -3
        UpperGeodesic.xmax= 3
        UpperGeodesic.ymax= 3

    ## plot commands
    @classmethod
    def set_plot_limits(cls):
        highest = max(abs(i)
                      for i in [cls.ymin, cls.ymax, cls.xmax, cls.xmin])
        cls.ax.axis([-highest, highest, 0, highest])

    @classmethod
    def plot_all(cls):
        if UpperGeodesic.ymax <= UpperGeodesic.ymin:
            UpperGeodesic.ymax = UpperGeodesic.ymin + 1 # else nothing to plot
        super().plot_all()
