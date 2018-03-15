#!/usr/bin/env python3

# class file geodesic.py

# started as a script to visualize what happens to hyperbolic plane
# if different isometries act on it

import geodesic as gd
import numpy as np
import numpy.linalg as lina
import matplotlib.pyplot as plt

# upper half space as the basic model

class DiscGeodesic(gd.Geodesic):
    """Geodesic line in disc model

    takes endpoints on boundary as arguments
    stores x and y data as points in x and y
    """
    
    def __init__(self, a, b, color="b", label = ''):
        """initialize Geodesic by endpoints

        parametrised by angle
        """

        super().__init__(a, b, color, label = '')

    def sort_se(self):
        """sort start and end
        
        self.start \in [0, 2 pi)
        end comes less than or equal to pi after start"""
        self.start = np.mod(self.start, 2*np.pi)
        self.end = np.mod(self.end, 2*np.pi)

        d = self.end - self.start
        if d > np.pi:
            d = 2*np.pi - d
            self.start = self.end
            self.end = self.start + d
        elif d < -np.pi:
            self.end += 2*np.pi
        elif d < 0:
            d = -d
            self.start = self.end
            self.end = self.start + d

    def get_data(self):
        eps = np.pi/100000
        if abs(self.end - self.start - np.pi) < eps:
            xs = np.linspace(np.cos(self.start), np.cos(self.end), self.res)
            ys = np.linspace(np.sin(self.start), np.sin(self.end), self.res)
        else:
            a = np.pi/2 - self.start
            b = a + (np.pi - self.end + self.start)
            cangle = (self.start + self.end)/2
            r = np.tan((self.end - self.start)/2)
            R = np.sqrt(1+r**2)
            xs = R*np.cos(cangle) + r*np.cos(np.linspace(a, b, self.res))
            ys = R*np.sin(cangle) + -r*np.sin(np.linspace(a, b, self.res))
        return(xs, ys)

    ####### TODO TODO ab hier noch nix gemacht!

    ## apply transformations to ONE geodesic
    def translate_one_geod(self, dx):
        a = self.from_disc_to_upper(self.start)
        b = self.from_disc_to_upper(self.end)
        if a != self.inf:
            a += dx
        if b != self.inf:
            b += dx
        self.start = self.from_upper_to_disc(a)
        self.end = self.from_upper_to_disc(b)
        self.sort_se()

    def translate_one_at_zero(self, dx):
        """change names to be more ... correct"""
        self.translate_one_geod(dx)

    def rotate_one_geod(self, phi):
        """rotates the geodesic around origin"""
        self.start += phi
        self.end += phi
        if self.start > 2*np.pi:
            self.sort_se()

    def hyperbolic_translate_one(self, dmult=1.001):
        a = self.from_disc_to_upper(self.start)
        b = self.from_disc_to_upper(self.end)
        if a != self.inf:
            a *= dmult
        if b != self.inf:
            b *= dmult
        self.start = self.from_upper_to_disc(a)
        self.end = self.from_upper_to_disc(b)
        self.sort_se()


    # tesselate hyperbolic space
    @classmethod
    def tesselate(cls, depth=5):
        """Tesselates the Durham hyperbolic plane"""
        start_angles = 2/3*np.pi*np.array([0, 1, 2])
        #start_angles = [0, np.pi/2, np.pi+0.1]
        new_triangles = []
        for i in range(len(start_angles)):
            a = start_angles[i]
            b = start_angles[(i+1)%3]
            c = start_angles[(i+2)%3]
            x = cls.mirror(a, b, c)
            cls(a, b, 'r')
            cls(a, x, 'orange')
            cls(x, b, 'orange')
            new_triangles.append( (a, b, x) )

        use_colors = ['yellow', 'g', 'b', 'indigo', 'violet']
        for i in range(depth):
            new_new_tr = []
            for t in new_triangles:
                a, b, x = t
                # order matters! (?) - ne...
                #n1 = cls.mirror(x, a, b)
                # n1 has some issues...
                n1 = cls.mirror(x, a, b, True)
                n2 = cls.mirror(b, x, a)
                cls(a, n1, use_colors[i])
                cls(x, n1, use_colors[i])
                cls(b, n2, use_colors[i])
                cls(x, n2, use_colors[i])
                new_new_tr.append( (a, x, n1) )
                new_new_tr.append( (x, b, n2) )
            new_triangles = new_new_tr


    ## plot commands
    @classmethod
    def plot_all(cls):
        super().plot_all()
        t = np.linspace(0, 2*np.pi, cls.res)
        xs = np.cos(t)
        ys = np.sin(t)
        #cls.ax.plot(xs, ys, color= 'k')
        cls.ax.plot(xs, ys, linewidth=3, color= 'k')

    # staticmethods
    @staticmethod
    def mirror(alpha, beta, x, has_issue = False):
        """mirrors x at the geodesic from a to b

        (everything given in radians)"""
        # TODO scheint falsch zu sein :(
        r = 1/np.cos((beta-alpha)/2)
        a = r*np.cos((alpha+beta)/2) - np.cos(x)
        b = r*np.sin((alpha+beta)/2) - np.sin(x)

        numerator = (a*np.cos(x) + b*np.sin(x))
        denom = (( a**2 + b**2 )**(1/2))
        t = numerator/denom
        if t < -1:
            t = -1 # somehow t was slightly less than -1 once...
        if has_issue:
            return np.mod(x + 2*np.arcsin(t), 2*np.pi)
        else:
            return np.mod(x - 2*np.arcsin(t), 2*np.pi)
