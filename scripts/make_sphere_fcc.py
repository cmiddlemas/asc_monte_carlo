#!/usr/bin/env python3
"""
make_sphere_fcc.py

Makes a fcc sphere packing at [density_mult]*phi_fcc with
[n_repeat]**3 particles. Sphere radius will be 1, and the output
will be in the form of asc_monte_carlo.

Usage: ./make_sphere_fcc.py [density_mult] [n_repeat] > outfile

Options:
    -h, --help: Display this message
"""

import sys
import numpy as np

if __name__ == "__main__":
    if "-h" in sys.argv or "--help" in sys.argv:
        print(__doc__)
        sys.exit()
    density_mult = float(sys.argv[1])
    n_repeat = int(sys.argv[2])

    scale = n_repeat*np.power(1.0/density_mult, 1.0/3.0)
    # https://mathworld.wolfram.com/CubicClosePacking.html
    # http://physics.bu.edu/~okctsui/PY543/3_notes_Crystals_2013.pdf
    u1 = scale*np.sqrt(2.0)*np.array([1.0, 1.0, 0.0])
    u2 = scale*np.sqrt(2.0)*np.array([1.0, 0.0, 1.0])
    u3 = scale*np.sqrt(2.0)*np.array([0.0, 1.0, 1.0])
    
    print("3 1 Sphere")
    print("{} {} {} {} {} {} {} {} {}".format(u1[0], u1[1], u1[2],
                                              u2[0], u2[1], u2[2],
                                              u3[0], u3[1], u3[2]))

    for i in range(n_repeat):
        for j in range(n_repeat):
            for k in range(n_repeat):
                v = np.array([float(i), float(j), float(k)])
                v /= float(n_repeat)
                print("{} {} {} 1.0".format(v[0], v[1], v[2]))
