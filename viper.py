#!/usr/bin/env python3

# starts as a script to visualize what happens to hyperbolic plane
# if different isometries act on it

# upper half space as the basic model

import uppergeodesic as ugd
import discgeodesic as dgd
import itertools as it
import numpy as np
import numpy.linalg as lina
import matplotlib.pyplot as plt
import matplotlib.animation as animation

import argparse

# argument handling
parser = argparse.ArgumentParser()
this_model = parser.add_mutually_exclusive_group()
this_model.add_argument('-d', '--disc', action='store_const', dest='model',
                        const='disc')
this_model.add_argument('-u', '--upper-half-space', action='store_const',
                        dest='model', const='uhs')

do_this = parser.add_mutually_exclusive_group(required=True)
do_this.add_argument('-l', '--hyperbolic', dest='ani_type',
                     action='store_const', const='hyperbolic')
do_this.add_argument('-p', '--parabolic', dest='ani_type',
                     action='store_const', const='parabolic')
do_this.add_argument('-e', '--elliptic', dest='ani_type',
                     action='store_const', const='elliptic')
do_this.add_argument('-t', '--tesselate', action='store_true')
do_this.add_argument('-s', '--subgroup', action='store_true')

parser.add_argument('-o', '--output', default=None, help='save the animation as an mp4 file')
parser.add_argument('-f', '--frames', default=201, help='number of frames')

args = parser.parse_args()

which_geodesic = {'disc': dgd.DiscGeodesic, 'uhs': ugd.UpperGeodesic}
Geodesic = which_geodesic.get(args.model, ugd.UpperGeodesic)


# different default functions

def main_animation(isom_type='hyperbolic'):
    # animates an isometry in hyperbolic space
    Geodesic.tesselate(4)
    this_animation = {
        'hyperbolic': (Geodesic.animate_hyperbolic_translation, 1.01),
        'parabolic': (Geodesic.animate_translation, 1/(args.frames-1)),
        'elliptic': (Geodesic.animate, np.pi / (args.frames-1))
    }
    ani_type = this_animation[isom_type]

    Geodesic.plot_all()
    ani = animation.FuncAnimation(Geodesic.fig,
                            ani_type[0],
                            frames=args.frames, fargs=(ani_type[1],), repeat=True,
                            interval=50)

    if args.output is None:
        plt.show()
    else:
        ani.save(args.output + ".mp4")


def main_tesselation():
    Geodesic.tesselate()
    Geodesic.plot_all()
    Geodesic.clean_plot()
    Geodesic.ax.axis([-1.01, 1.01, -1.01, 1.01])
    Geodesic.fig.set_size_inches(10, 10)
    plt.savefig("tesselation.png", bbox_inches='tight', pad_inches=0,
                transparent=True)
    plt.show()


def subgroup_tesselation():
    # generators of the group
    # T = np.matrix([[1,1],[0,1]])
    # S = np.matrix([[0,1],[-1,0]])

    # Kongruenzuntergruppe (wiki seite lemma von selberg)
    # T2 = T**3
    # STS = S*T2*S
    # T2i = lina.inv(T2)
    # STSi = lina.inv(STS)
    # print(T2)
    # print(STS)

    # generators of the subgroup
    gens = [np.matrix([[1, 2], [0, 1]]), np.matrix([[1, 0], [2, 1]]),
            np.matrix([[1, 2], [2, 5]]), np.matrix([[1, -2], [-2, 5]]),
            np.matrix([[1, -2], [2, -3]]), np.matrix([[1, 2], [-2, -3]]),
            np.matrix([[5, 2], [2, 1]]), np.matrix([[5, -2], [-2, 1]]),
            np.matrix([[-3, -2], [2, 1]]), np.matrix([[-3, 2], [-2, 1]]),
            np.matrix([[3, 2], [4, 3]]), np.matrix([[3, 4], [2, 3]])]
    inv_gens = [lina.inv(A) for A in gens]

    g0 = Geodesic(1, Geodesic.inf)
    g1 = Geodesic(-1, Geodesic.inf)
    g2 = Geodesic(-1, 0)
    g3 = Geodesic(1, 0)
    start_geods = [g0, g1, g2, g3]
    tmp = []
    tmpi = []
    new_geods = [g.new_from_matrix(A) for A in gens for g in start_geods]
    new_from_inverse = [g.new_from_matrix(A)
                        for A in inv_gens for g in start_geods]
    for i in range(1):
        new_geods += new_from_inverse
        for geod in new_geods:
            tmp += [geod.new_from_matrix(A) for A in gens]
            tmpi += [geod.new_from_matrix(A) for A in inv_gens]
        new_geods = tmp
        new_from_inverse = tmpi
        tmp = []
        tmpi = []

    Geodesic.xmin = -2
    Geodesic.xmax = 2
    Geodesic.ymax = 2
    Geodesic.plot_all()
    Geodesic.fig.set_size_inches(40, 20)
    plt.axis([-3.5, 3.5, 0, 2])
    plt.savefig("subgroup.pdf")
    plt.show()


def group():
    # print matrices in group generated by stuff
    t = np.matrix([[1, 1], [0, 1]])
    s = np.matrix([[0, 1], [-1, 0]])

    t2 = t ** 2
    sts = s * t * s
    t2i = lina.inv(t2)
    stsi = lina.inv(sts)
    new_matrices = [t2, sts, t2i, stsi]
    generators = [t2, sts, t2i, stsi]
    tmp = []
    for i in range(3):
        for m in new_matrices:
            for g in generators:
                tmp.append(m * g)
                print(tmp[-1])
        new_matrices = tmp
        tmp = []


def disc():
    dgd.DiscGeodesic(np.pi / 4, np.pi, 'r')
    dgd.DiscGeodesic(np.pi / 5, np.pi, 'g')
    dgd.DiscGeodesic(np.pi / 2, np.pi * 3 / 2, 'k')
    dgd.DiscGeodesic.plot_all()
    plt.show()


if args.tesselate:
    main_tesselation()
elif args.subgroup:
    subgroup_tesselation()
else:
    main_animation(args.ani_type)
