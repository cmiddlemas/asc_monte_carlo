#!/usr/bin/env python3
"""
Script to animate a series of saved sphere configs

Based on example from
https://matplotlib.org/3.1.1/gallery/animation/rain.html#sphx-glr-gallery-animation-rain-py

Usage: ./animate.py [files]

Options:
    --help: print this message
"""

import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.patches import Circle
from matplotlib.lines import Line2D 

def update(frame):
    print(frame)
    try:
        f = open(sys.argv[frame+1])
        f.readline()
        s = f.readline()
        unit = np.fromstring(s, sep = ' ')
        config = np.atleast_2d(np.loadtxt(f))
        ax.clear()
        ax.set_xlim(-3.0 -config_first[0,2]*1.1, unit_first[0]+3.0+config_first[0,2]*1.1)
        ax.set_ylim(-3.0 -config_first[0,2]*1.1, unit_first[3]+3.0+config_first[0,2]*1.1)
        # https://stackoverflow.com/questions/36470343/how-to-draw-a-line-with-matplotlib
        #print(unit)
        for i in range(0,len(config)):
            circle = Circle((config[i,0], config[i,1]), config[i,2])
            ax.add_artist(circle)
        line_w = 2.0
        ax.add_line(Line2D([0.0,unit[0]],[0.0,unit[1]], color="black",linewidth=line_w))
        ax.add_line(Line2D([0.0,unit[2]],[0.0,unit[3]], color="black", linewidth=line_w))
        ax.add_line(Line2D([unit[0],unit[0]+unit[2]],[unit[1],unit[1]+unit[3]], color="black", linewidth=line_w))
        ax.add_line(Line2D([unit[2],unit[0]+unit[2]],[unit[3],unit[1]+unit[3]], color="black", linewidth=line_w))
    except:
        return
    
    return


if __name__ == "__main__":
    if "--help" in sys.argv:
        print(__doc__)
        sys.exit(0)

    # https://stackoverflow.com/questions/332289/how-do-you-change-the-size-of-figures-drawn-with-matplotlib
    fig = plt.figure(figsize=(5,5))
    ax = fig.add_axes([0, 0, 1, 1], frameon=True)
    ax.axis('equal')

    f = open(sys.argv[1])
    f.readline()
    s = f.readline()
    unit_first = np.fromstring(s, sep = ' ')
    config_first = np.atleast_2d(np.loadtxt(f))


    # https://stackoverflow.com/questions/33181699/matplotlib-animation-not-saving-properly
    animation = FuncAnimation(fig, update, interval=100, frames=len(sys.argv)-1, repeat=False)
    # https://towardsdatascience.com/animations-with-matplotlib-d96375c5442c
    animation.save('animate.gif', writer='imagemagick')
