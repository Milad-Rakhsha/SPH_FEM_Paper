import csv,os,sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import pandas as pd
from decimal import Decimal
from collections import OrderedDict
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'


path = str(sys.argv[1])

Yang_x=path+"/YangEtAl10X.txt"
Yang_y=path+"/YangEtAl10Y.txt"
Antoci_x=path+"/Antoci_x.txt"
Antoci_y=path+"/Antoci_y.txt"
Ex_x=path+"/Ex10X.txt"
Ex_y=path+"/Ex10Y.txt"
FEM="591793.csv"
FEM1="591794.csv"
FEM2="591795.csv"
FEM3="591796.csv"
FEM4="591797.csv"
FEM5="591798.csv"
FEM6="591799.csv"
FEM7="591800.csv"
FEM8="591801.csv"

SPH="591190.txt"
SPH1="591191.txt"
SPH2="591192.txt"
SPH3="591193.txt"

files = OrderedDict([#
    (Yang_y ,{"axis":"y",   "shift":-2.011,"dT": 0.02, "lineStyle": 'b--', "markerEvery": 1, 'markersize':5, 'label': 'Yang et al.'}),
    (Antoci_y ,{"axis":"y",   "shift":-2.011,"dT": 0.02, "lineStyle": 'r--', "markerEvery": 1, 'markersize':5, 'label': 'Antoci et al.'}),
    (Ex_y ,{"axis":"y",   "shift":-2.011,"dT": 0.02, "lineStyle": 'ko', "markerEvery": 1, 'markersize':5, 'label': 'Experimental'}),
    # (FEM,  {"axis": "y", "shift": +0.65, "dT": 0.05, "lineStyle": 'g-', "markerEvery": 1, 'markersize':5, 'label': 'FEM'}),
    (FEM1,  {"axis": "y", "shift": +0.65, "dT": 0.05, "lineStyle": '-', "markerEvery": 1, 'markersize':5, 'label': 'FEM'}),
    # (FEM2,  {"axis": "y", "shift": +0.65, "dT": 0.05, "lineStyle": '-', "markerEvery": 1, 'markersize':5, 'label': 'FEM2'}),
    # (FEM3,  {"axis": "y", "shift": +0.65, "dT": 0.05, "lineStyle": '-', "markerEvery": 1, 'markersize':5, 'label': 'FEM3'}),
    # (FEM4,  {"axis": "y", "shift": +0.65, "dT": 0.05, "lineStyle": '-', "markerEvery": 1, 'markersize':5, 'label': 'FEM4'}),
    # (FEM5,  {"axis": "y", "shift": +0.65, "dT": 0.05, "lineStyle": '-', "markerEvery": 1, 'markersize':5, 'label': 'FEM5'}),
    # (FEM6,  {"axis": "y", "shift": +0.65, "dT": 0.05, "lineStyle": '-', "markerEvery": 1, 'markersize':5, 'label': 'FEM6'}),
    # (FEM7,  {"axis": "y", "shift": +0.65, "dT": 0.05, "lineStyle": '-', "markerEvery": 1, 'markersize':5, 'label': 'FEM7'}),

    (SPH,  {"axis":  "z",  "shift": +0.65, "dT": 0.05, "lineStyle": '-.', "markerEvery": 1, 'markersize':5, 'label': 'SPH'}),
    # (SPH1,  {"axis": "z", "shift": +0.65, "dT": 0.05, "lineStyle": '-.', "markerEvery": 1, 'markersize':5, 'label': 'SPH1'}),
    # (SPH2,  {"axis": "z", "shift": +0.65, "dT": 0.05, "lineStyle": '-.', "markerEvery": 1, 'markersize':5, 'label': 'SPH2'}),
    # (SPH3,  {"axis": "z", "shift": +0.65, "dT": 0.05, "lineStyle": '-.', "markerEvery": 1, 'markersize':5, 'label': 'SPH3'}),

])

major_ticks = np.arange(0, 101, 20)
minor_ticks = np.arange(0, 101, 5)


fig = plt.figure(num=None,figsize=(12, 8),  dpi=300, facecolor='w', edgecolor='k')
ax1 = fig.add_subplot(212)
ax1.set_ylabel('y($mm$)')
ax1.grid(color='k', linestyle='-', linewidth=0.2)
ax1.autoscale(enable=True, axis='x', tight=True)

for thisFile in files:
    data = pd.read_csv(thisFile)
    axis_=files[thisFile]["axis"]
    time=data["Time"]
    x =(data[axis_]-data[axis_][0])*1000
    ax1.plot(time, x,
              files[thisFile]["lineStyle"],
              linewidth=1, markersize=files[thisFile]["markersize"], markevery=files[thisFile]["markerEvery"],
              label=files[thisFile]['label'])



plt.legend( fancybox=True, shadow=True, ncol=1)
ax1.set_xticks(np.linspace(0, 0.4, 9))
ax1.set_yticks(np.linspace(0, 20, 11))
ax1.grid(which='both', linestyle='--', linewidth=0.5)
ax1.set_xlim(0, 0.4)
ax1.set_ylim(0, 18)
ax1.set_yticks(np.linspace(0, 20, 21), minor=True)
ax1.grid(which='minor', alpha=0.2)
ax1.grid(which='major', alpha=0.3)



files = OrderedDict([#
    (Yang_x,{"axis": "x", "shift": +0.65, "dT": 0.05, "lineStyle": 'b--', "markerEvery": 1, 'markersize':5, 'label': 'Yang et al.'}),
    (Antoci_x,{"axis": "x", "shift": +0.65, "dT": 0.05, "lineStyle": 'r--', "markerEvery": 1, 'markersize':5, 'label': 'Antoci et al.'}),
    (Ex_x,{"axis": "x", "shift": +0.65, "dT": 0.05, "lineStyle": 'ko', "markerEvery": 1, 'markersize':5, 'label': 'Experimental'}),
    # (FEM,{"axis": "x", "shift": +0.65, "dT": 0.05, "lineStyle": '-', "markerEvery": 1, 'markersize':5, 'label': 'FEM'}),
    (FEM1,{"axis": "x", "shift": +0.65, "dT": 0.05, "lineStyle": '-', "markerEvery": 1, 'markersize':5, 'label': 'FEM'}),
    # (FEM2,{"axis": "x", "shift": +0.65, "dT": 0.05, "lineStyle": '-', "markerEvery": 1, 'markersize':5, 'label': 'FEM2'}),
    # (FEM3,{"axis": "x", "shift": +0.65, "dT": 0.05, "lineStyle": '-', "markerEvery": 1, 'markersize':5, 'label': 'FEM3'}),
    # (FEM4,{"axis": "x", "shift": +0.65, "dT": 0.05, "lineStyle": '-', "markerEvery": 1, 'markersize':5, 'label': 'FEM4'}),
    # (FEM5,{"axis": "x", "shift": +0.65, "dT": 0.05, "lineStyle": '-', "markerEvery": 1, 'markersize':5, 'label': 'FEM5'}),
    # (FEM6,{"axis": "x", "shift": +0.65, "dT": 0.05, "lineStyle": '-', "markerEvery": 1, 'markersize':5, 'label': 'FEM6'}),
    # (FEM7,{"axis": "x", "shift": +0.65, "dT": 0.05, "lineStyle": '-', "markerEvery": 1, 'markersize':5, 'label': 'FEM7'}),
    # (FEM8,{"axis": "x", "shift": +0.65, "dT": 0.05, "lineStyle": '-', "markerEvery": 1, 'markersize':5, 'label': 'FEM8'}),

    (SPH,   {"axis": "x", "shift": +0.65, "dT": 0.05, "lineStyle": '-.', "markerEvery": 1, 'markersize':5, 'label': 'SPH'}),
    # (SPH1,  {"axis": "x", "shift": +0.65, "dT": 0.05, "lineStyle": '-.', "markerEvery": 1, 'markersize':5, 'label': 'SPH1'}),
    # (SPH2,  {"axis": "x", "shift": +0.65, "dT": 0.05, "lineStyle": '-.', "markerEvery": 1, 'markersize':5, 'label': 'SPH2'}),
    # (SPH3,  {"axis": "x", "shift": +0.65, "dT": 0.05, "lineStyle": '-.', "markerEvery": 1, 'markersize':5, 'label': 'SPH3'}),
])

ax2 = fig.add_subplot(211)
ax2.grid(color='k', linestyle='-', linewidth=0.2)
ax2.set_ylabel('x($mm$)')
ax2.set_xlabel("time(s)")
ax2.set_title("FSI Validation")

ax2.autoscale(enable=True, axis='x', tight=True)
for thisFile in files:
    data = pd.read_csv(thisFile)
    axis_=files[thisFile]["axis"]
    time=data["Time"]
    x =(data[axis_]-data[axis_][0])*1000
    ax2.plot(time, x,
              files[thisFile]["lineStyle"],
              linewidth=1, markersize=files[thisFile]["markersize"], markevery=files[thisFile]["markerEvery"],
              label=files[thisFile]['label'])

leg = ax1.legend()
ax2.set_xticks(np.linspace(0, 0.4, 9))
ax2.set_yticks(np.linspace(0, 50, 11))
ax2.grid(which='both', linestyle='--', linewidth=0.5)
ax2.set_xlim(0, 0.4)
ax2.set_ylim(0, 20)
ax2.set_yticks(np.linspace(0, 50, 21), minor=True)
ax2.grid(which='minor', alpha=0.2)
ax2.grid(which='major', alpha=0.3)
plt.savefig('Fig.png', bbox_inches='tight')
plt.show()
