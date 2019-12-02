#!/usr/bin/env python3
#
# Purpose:  read torque time series from ascii file, compute average, fft, plot
#           everything, and print difference between inner and outer torque
# Usage:    python torque.py 
# Author:   Daniel Feldmann
# Date:     11th June 2018
# Modified: 14th March 2019
#
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['text.usetex'] = True
mpl.rcParams['text.latex.preamble'] = [
    r"\usepackage[utf8x]{inputenc}",
    r"\usepackage[T1]{fontenc}",
    r"\usepackage[detect-all]{siunitx}",
    r'\usepackage{amsmath, amstext, amssymb}',
    r'\usepackage{xfrac}',
    r'\usepackage{lmodern, palatino, eulervm}']
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams.update({'font.size': 10})

# read time series data from ascii file
fnam = 'torque'
print('Read time series data from file', fnam)
f = open(fnam, 'r')
#t, ti, to = np.genfromtxt(f, unpack=True)
t, ti, to, retaui, retauo = np.genfromtxt(f, unpack=True)
f.close()
print('Data from t =', t[0], 'to t =', t[-1])

# set time window for averaging and fft
ta0 = t[0]   # default start time
ta1 = t[-1]  # default end time 
ta0 =  10.0  # set later start time
#ta1 = 20.0 # set end time
nta = np.where((t >= ta0) & (t<= ta1)) # index range for averaging
ta  = t[nta] 
ta0 = np.min(ta) # reset time window to actual bounds
ta1 = np.max(ta)  
print ('Temporal averaging from ta0 =', ta0, ' to ta1 =', ta1)

# average torque
tim = np.mean(ti[nta])
tom = np.mean(to[nta])
print('Averaged inner torque <ti>_ta =', tim)
print('Averaged outer torque <to>_ta =', tom)
print('Absolute torque difference  o-i    = ',  tom-tim)
print('Relative torque difference (o-i)/o = ', (tom-tim)/tim)
timp = np.ones(2) * tim # saving memory plotting mean over relevant range
tomp = np.ones(2) * tom
tap  = np.zeros(2)
tap[0] = ta[0]
tap[1] = ta[-1]

# torque fluctuation
tif = ti[nta]-tim
tof = to[nta]-tom
del nta # free memory

# create figure suitable for A4 format
def mm2inch(*tupl):
 inch = 25.4
 if isinstance(tupl[0], tuple):
  return tuple(i/inch for i in tupl[0])
 else:
  return tuple(i/inch for i in tupl)
#fig = plt.figure(num=None, figsize=mm2inch(210.0,297.0), dpi=600)
fig = plt.figure(num=None, figsize=mm2inch(160.0,140.0), dpi=150)

# plot inner and outer torque over time
ax1 = plt.subplot2grid((3, 1), (0, 0), rowspan=1, colspan=1)
#ax1.set_title(r"Temporal evolution of torque, $Ra=3550$, \texttt{tc0076}")
plt.setp(ax1.get_xticklabels(), visible=False)
ax1.set_ylabel(r"$N\!u_{\omega}$")
ax1.plot(t,   ti,   color='blue', linestyle='-',  label=r"Inner")
ax1.plot(tap, timp, color='black', linestyle='--', label=r"$\langle N\!u_{\omega,i} \rangle_{t_a}$")
ax1.plot(t,   to,   color='red',  linestyle='-',  label=r"Outer")
ax1.plot(tap, tomp, color='black', linestyle=':', label=r"$\langle N\!u_{\omega,o} \rangle_{t_a}$")
ax1.legend(loc='lower right', ncol=4)
del t, ti, to # free memory

# plot inner torque fluctuation
ax2 = plt.subplot2grid((3, 1), (1, 0), rowspan=1, colspan=1, sharex=ax1)
ax2.set_ylabel(r"$N\!u_{\omega}$")
plt.setp(ax2.get_xticklabels(), visible=False)
ax2.plot(ta, tif, color='blue', linestyle='-', label=r"$N\!u_{\omega,i}^{\prime}=N\!u_{\omega,i}-\langle N\!u_{\omega,i} \rangle_{t_a}$")
ax2.legend(loc='best', ncol=1)

# plot outer torque fluctuation
ax3 = plt.subplot2grid((3, 1), (2, 0), rowspan=1, colspan=1, sharex=ax1)
ax3.set_ylabel(r"$N\!u_{\omega}$")
ax3.set_xlabel(r"$t$ in $\sfrac{d}{\nu^{2}}$")
ax3.plot(ta, tof, color='red', linestyle='-', label=r"$N\!u_{\omega,o}^{\prime}=N\!u_{\omega,o}-\langle N\!u_{\omega,o} \rangle_{t_a}$")
ax3.legend(loc='best', ncol=1)
del ta, tif, tof # free memory

# io mode
pdf = 1
if pdf != 1:
 # interactive output
 plt.tight_layout()
 plt.show()
else:
 # pdf file output
 #outf = "plotLongSection{:04d}.pdf".format(n)
 fig.tight_layout()
 outf = "torque.pdf"
 plt.savefig(outf)
 print("Written file", outf)
 fig.clf()
