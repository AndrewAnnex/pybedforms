import argparse
import numpy as np
import configparser
from .dune_topo import DuneTopo


def main(color, k, n, f, r, t, w, erosion, d, *infiles):
    c = configparser.ConfigParser()
    for file in infiles:
        with open(file) as src:
            c.read_file(src)
        config = c.defaults()
        make_movie(color, k, n, f, r, t, w, erosion, d, config)
    pass

def make_movie(color, k, n, f, r, t, w, erosion, d, config):
    m = k + n
    q = m + f
    j = r / f
    s = q + t
    v = s + w
    b = -erosion * w
    p = v + d
    if color in ('color','c','Color','COLOR'):
        Total_Frames_Per_Movie = p
    #
    dT = np.arange(0, k)
    dT[k + 1: m] = k
    dT[m + 1: p] = k

    dTrend = np.arange(0, m)
    dTrend[m + 1: q] = np.arange(j, j, r)
    dTrend[q + 1: s] = r
    dTrend[s + 1: p] = r

    dZHO = np.arange(0,s)
    dZHO[s + 1: v] = np.arange(-erosion, -erosion, b)
    dZHO[v + 1: p] = b
    NumberOfFrames = len(dT)

    # make the dune topo
    dune = DuneTopo(**config)

    pass


if __name__ == '__main__':
    p = argparse.ArgumentParser(description='New Dune Tool')
    """
    How many files do you want to run? 
    Name of input parameter file?
    Name of input parameter file?
    Enter “color” for movies/tiffs or “bw” for single image post-scripts :  color
    How many frames showing deposition?  : 100
    How many pause frames after deposition?  : 5
    How many frames showing rotation?  : 90
    How many rotational degrees?  : 45
    How many pause frames after rotation?  : 5
    How many frames showing erosion?  : 100
    How much erosion between frames?  : 0.025
    How many pause frames after erosion? 
    """
    p.add_argument('color', default='bw', type=str, dest='color')
    p.add_argument('k', type=int, help='How many frames showing deposition?')
    p.add_argument('n', type=int, help='How many pause frames after deposition?')
    p.add_argument('f', type=int, help='How many frames showing rotation?')
    p.add_argument('r', type=int, help='How much total rotation in degrees?')
    p.add_argument('t', type=int, help='How many pause frames after rotation?')
    p.add_argument('w', type=int, help='How many frames showing erosion?')
    p.add_argument('erosion', type=float, help='How much erosion between frames? (0.025 suggested)', default=0.025)
    p.add_argument('d', type=int, help='How many pause frames at the end?')
    p.add_argument('infiles', nargs='+')
    args = p.parse_args()
    main(*args)