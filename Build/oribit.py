#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np
from scipy.io import FortranFile
f = FortranFile('sys.n_body', 'r','>u4')


n_bodies=f.read_ints('>u4')
calc_len=f.read_reals('>f8')
labels=f.read_record('a15')
time=0

prog=[]
t=[]
while time<calc_len:

    try:
        time=f.read_reals('>f8')
        print(time)
        pos=f.read_reals('>f8').reshape((int(n_bodies),3),order="F")
        vel=f.read_reals('>f8').reshape((int(n_bodies),3),order="F")
        t.append(time)
        prog.append(pos)
    except:
        break




for i in range(0,n_bodies[0]):
    x=[]
    y=[]
    for dt in range(len(prog)):

        x.append(prog[dt][i,0])
        y.append(prog[dt][i,1])

    plt.plot(x,y,label=labels[i])
if n_bodies<10:
    plt.legend()
plt.show()
