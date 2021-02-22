VIPer - Visualiser of Isometries in the Plane of hypERbolic geometry
====================================================================

DESCRIPTION
-----------

VIPer can draw tesselations of the hyperbolic plane in the disc or upper half
space model. It either just draws the tesselation or animates certain elliptic,
parabolic or loxodromic isometries. It can also draw the tesselation for the
principal congruence subgroup of level 2 of the modular group.

BASIC USAGE
-----------

```
usage: viper.py [-h] [-d | -u] (-l | -p | -e | -t | -s)

optional arguments:
  -h, --help            show this help message and exit
  -d, --disc
  -u, --upper-half-space
  -l, --hyperbolic
  -p, --parabolic
  -e, --elliptic
  -t, --tesselate
  -s, --subgroup
```

The options `-d` and `-u` choose the Poincaré disk or the upper half space model, respectively.
The options `--hyperbolic`, `--parabolic`, and `--elliptic` animate an appropriate isometry of the hyperbolic plane.
The options `--tesselate` and `--subgroup` do not show an animation but plot tesselations of the hyperbolic plane with a fundamental domain of the modular group or of a torsion free subgroup.

BONUS
-----

The files [pentagon.py](pentagon.py) and [hexagon.py](hexagon.py) implement a two-dimensional version of the algorithm explained in [this paper](https://doi.org/10.1007/s10711-018-0357-y).

TODO/IDEAS
----------

 - implement different tilings
 - implement more flexibility for the isometries
 - draw other subgroup tesselations
 - draw the dirichlet domain for a group given by generators (branch dirichlet)
