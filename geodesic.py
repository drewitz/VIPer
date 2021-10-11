#!/usr/bin/env python3

# class file geodesic.py

# started as a script to visualize what happens to hyperbolic plane
# if different isometries act on it

import numpy as np
import numpy.linalg as lina
import matplotlib.pyplot as plt

# upper half space as the basic model

class Geodesic:
    """Geodesic line in hyperbolic space

    takes endpoints on boundary as arguments
    stores x and y data as points in x and y
    string "inf" is point at infinity, i.e. y=inf
    """
    
    res = 1001 # resolution to plot later
    all_geods = [] # no geodesics yet
    plotted_geods = [] # plotted geodesics, gets reset anyway.
    inf = "inf"
    legend = False

    def __init__(self, a, b, color="b", label=''):
        """initialize Geodesic by endpoints

        a, b - x values of the endpoints or "inf" if infinity
        res is resolution
        """
        assert a != b, "Geodesic from point to itself is empty"

        self.color = color
        self.label = label
        if label != '':
            self.__class__.legend = True
        self.start = a
        self.end = b
        self.sort_se()
        # if one is inf, then 'self.start'
        # else start < end
        # hence -> unique order

        new = True
        for g in self.__class__.all_geods:
            if self.start == g.start and self.end == g.end:
                new = False
        if new:
            self.__class__.all_geods.append(self)

    def __comp__(self, other):
        # doesn't seem to work in __init__ :(
        if not isinstance(other, self.__class__):
            return(NotImplemented)
        if self.start == other.start and self.end == other.end:
            return(True)
        else:
            return(False)

    def __repr__(self):
        return f"Geodesic{self}"

    def __str__(self):
        return f"({self.start}, {self.end})"

    def sort_se(self):
        """sort start and end - depends on model"""
        raise NotImplementedError

    def get_data(self):
        raise NotImplementedError

    ## the next two functions create new geodesics from existing ones
    def new_geod(self, a, b, c, d):
        # only in upper half space
        raise NotImplementedError

    def new_from_matrix(self, M):
        # only in upper half space
        raise NotImplementedError

    ## apply transformations to ONE geodesic
    def apply_matrix(self, M):
        raise NotImplementedError

    def translate_one_geod(self, dx):
        raise NotImplementedError

    def translate_one_at_zero(self, dx):
        raise NotImplementedError

    def rotate_one_geod(self, phi):
        raise NotImplementedError

    def hyperbolic_translate_one(self, dmult=1.001):
        raise NotImplementedError

    # apply transformations to ALL geodesics
    @classmethod
    def hyperbolic_translate_all(cls, d=1.001):
        """translate hyperbolic plane along the Geodesic(-1,1)"""
        for g in cls.all_geods:
            g.hyperbolic_translate_one(d)

    @classmethod
    def rotate_all(cls, phi):
        for geod in cls.all_geods:
            geod.rotate_one_geod(phi)

    @classmethod
    def translate_all(cls, dx):
        for geod in cls.all_geods:
            geod.translate_one_at_zero(dx)

    # Animations
    @classmethod
    def animate_persistent(cls, dphi=np.pi/500):
        """animates with new draws instead of redraws"""
        cls.rotate_all(dphi)
        for geod in cls.all_geods:
            geod.plot_one_geod()

    @classmethod
    def animate_hyperbolic_translation(cls, dmult=1.001):
        cls.hyperbolic_translate_all(dmult)
        cls.update_plot()

    @classmethod
    def animate(cls, dphi=np.pi/180):
        assert len(cls.plotted_geods)==len(cls.all_geods),\
                "plotted {0} geodesic but I have {1}".format(
                        len(cls.plotted_geods), len(cls.all_geods))

        cls.rotate_all(dphi)
        cls.update_plot()

    @classmethod
    def animate_translation(cls, dx=0.01):
        assert len(cls.plotted_geods)==len(cls.all_geods),\
                "plotted {0} geodesic but I have {1}".format(
                        len(cls.plotted_geods), len(cls.all_geods))

        cls.translate_all(dx)
        cls.update_plot()

    # tesselate hyperbolic space
    @classmethod
    def tesselate(self, depth=10):
        raise NotImplementedError

    ## plot commands
    def plot_one_geod(self):
        """plots one geodesic

        used in plot_all and relies on cls.ax being defined!
        """
        xs, ys = self.get_data()
        line, = self.__class__.ax.plot(xs, ys, self.color, label=self.label)
        #line, = self.__class__.ax.plot(xs, ys, self.color, linewidth=3)
        self.__class__.plotted_geods.append(line)
        return(line)

    @classmethod
    def set_plot_limits(cls):
        raise NotImplementedError

    @classmethod
    def plot_all(cls):
        cls.fig, cls.ax = plt.subplots()
        cls.ax.set_aspect("equal")
        #cls.set_plot_limits() # TODO fuer uppergeod wieder reaktivieren
        cls.plotted_geods = []

        for geod in cls.all_geods:
            line = geod.plot_one_geod()
        if cls.legend:
            plt.legend(handles = cls.plotted_geods)

    @classmethod
    def update_plot(cls):
        for i, line in enumerate(cls.plotted_geods):
            xs, ys = cls.all_geods[i].get_data()
            line.set_data(xs, ys)

    @classmethod
    def clean_plot(cls):
        cls.ax.axis('off')

    # static methods as helper functions
    @staticmethod
    def inversion_on_unit_circle(x):
        """Inversion at unit circle

        maps upper half space to upper half space, swaps 0 and inf
        (x,0) -> (y,0)
        """
        y = 0
        if x != Geodesic.inf:
            eps = 10**(-6)
            if abs(x) < eps:
                y = Geodesic.inf
            else:
                y = 1/x
        return(y)

    @staticmethod
    def apply_moebius(a, b, c, d, x):
        """apply moebius transformation on x"""
        eps = 1/1000
        dest = Geodesic.inf
        if c != 0:
            dest = a/c

        if x != Geodesic.inf:
            if abs(c*x + d) < eps:
                dest = Geodesic.inf
            else:
                dest = (a*x + b)/(c*x + d)
        return(dest)

    @staticmethod
    def from_disc_to_upper(phi):
        """inversion from disc to upper half space - boundary

        (cos(phi), sin(phi)) -> (x,0)
        """
        while phi > np.pi:
            phi -= 2*np.pi
        while phi < -np.pi:
            phi += 2*np.pi
        # now -pi < phi <= pi
        eps = np.pi/100000 # this will be error prone...
        if abs(phi+np.pi/2) < eps:
            x = Geodesic.inf
        else:
            x = np.sqrt( (1-np.sin(phi))/(1+np.sin(phi)))
            if phi > np.pi/2 or phi < -np.pi/2:
                x *= -1
        return(x)

    @staticmethod
    def from_upper_to_disc(x):
        """inversion from upper half space to disc - boundary

        (x,0) -> (cos(phi), sin(phi))
        """
        phi = -np.pi/2
        if x != Geodesic.inf:
            phi = np.arcsin((1-x**2)/(1+x**2))
            if x < 0:
                # TODO check this
                phi = np.pi - phi
        return(phi)
