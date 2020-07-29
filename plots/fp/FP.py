from decimal import Decimal
import matplotlib
matplotlib.use('TkAgg')
# matplotlib.use('Agg')  # Must be before importing matplotlib.pyplot or pylab!
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from collections import OrderedDict
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
matplotlib.rcParams.update({'font.size': 18})

import scipy.signal as signal


def FEMFunction(vect):
    # First, design the Buterworth filter
    N = 5    # Filter order
    Wn = 0.5  # Cutoff frequency
    B, A = signal.butter(N, Wn, output='ba')
    return signal.filtfilt(B, A, vect)


fem_1=("590573.csv")
fem_2=("590567.csv")


isph_1=("1")
isph_1=("2")

files = OrderedDict([#
    (fem_1,{"method": "FEM", "F":"Fx", "t":"Time", "dt": 0.0015, "lineStyle": '--', "markerEvery": 1, 'markersize':5, 'label': 'FEM-1'}),
    (fem_2,{"method": "FEM", "F":"Fx", "t":"Time", "dt": 0.0015, "lineStyle": '--', "markerEvery": 1, 'markersize':5, 'label': 'FEM-2'}),
    # (fem_3,{"method": "FEM", "F":"Fy", "t":"t", "dt": 0.0015, "lineStyle": '--', "markerEvery": 1, 'markersize':5, 'label': 'FEM-3'}),
    # (fem_4,{"method": "FEM", "F":"Fy", "t":"t", "dt": 0.0015, "lineStyle": '--', "markerEvery": 1, 'markersize':5, 'label': 'FEM-4'}),
    # (fem_5,{"method": "FEM", "F":"Fy", "t":"t", "dt": 0.0015, "lineStyle": '--', "markerEvery": 1, 'markersize':5, 'label': 'FEM-5'}),
    # (fem_6,{"method": "FEM", "F":"Fy", "t":"t", "dt": 0.0015, "lineStyle": '--', "markerEvery": 1, 'markersize':5, 'label': 'FEM-6'}),
    # (fem_7,{"method": "FEM", "F":"Fy", "t":"t", "dt": 0.0015, "lineStyle": '--', "markerEvery": 1, 'markersize':5, 'label': 'FEM-7'}),
    # (fem_8,{"method": "FEM", "F":"Fy", "t":"t", "dt": 0.0015, "lineStyle": '--', "markerEvery": 1, 'markersize':5, 'label': 'FEM-8'}),
    # (fem_9,{"method": "FEM", "F":"Fy", "t":"t", "dt": 0.0015, "lineStyle": '--', "markerEvery": 1, 'markersize':5, 'label': 'FEM-9'}),
    # (fem_10,{"method": "FEM", "F":"Fy", "t":"t", "dt": 0.0015, "lineStyle": '--', "markerEvery": 1, 'markersize':5, 'label': 'FEM-10'}),
    # (fem_11,{"method": "FEM", "F":"Fy", "t":"t", "dt": 0.0015, "lineStyle": '--', "markerEvery": 1, 'markersize':5, 'label': 'FEM-11'}),

    # (fem_12,{"method": "FEM", "F":"Fy", "t":"t", "dt": 0.0015, "lineStyle": '--', "markerEvery": 1, 'markersize':5, 'label': 'FEM-12'}),
    # (fem_13,{"method": "FEM", "F":"Fy", "t":"t", "dt": 0.0015, "lineStyle": '--', "markerEvery": 1, 'markersize':5, 'label': 'FEM-13'}),
    # (fem_14,{"method": "FEM", "F":"Fy", "t":"t", "dt": 0.0015, "lineStyle": '--', "markerEvery": 1, 'markersize':5, 'label': 'FEM-14'}),
    # (fem_15,{"method": "FEM", "F":"Fy", "t":"t", "dt": 0.0015, "lineStyle": '--', "markerEvery": 1, 'markersize':5, 'label': 'FEM-15'}),
    # (fem_17,{"method": "FEM", "F":"Fy", "t":"t", "dt": 0.0015, "lineStyle": '--', "markerEvery": 1, 'markersize':5, 'label': 'FEM-17'}),
    # (fem_18,{"method": "FEM", "F":"Fy", "t":"t", "dt": 0.0015, "lineStyle": '--', "markerEvery": 1, 'markersize':5, 'label': 'FEM-18'}),
    # (fem_19,{"method": "FEM", "F":"Fy", "t":"t", "dt": 0.0015, "lineStyle": '--', "markerEvery": 1, 'markersize':5, 'label': 'FEM-19'}),
    # (fem_20,{"method": "FEM", "F":"Fy", "t":"t", "dt": 0.0015, "lineStyle": '--', "markerEvery": 1, 'markersize':5, 'label': 'FEM-20'}),
    # (fem_21,{"method": "FEM", "F":"Fy", "t":"t", "dt": 0.0015, "lineStyle": '--', "markerEvery": 1, 'markersize':5, 'label': 'FEM-21'}),
    # (fem_22,{"method": "FEM", "F":"Fy", "t":"t", "dt": 0.0015, "lineStyle": '--', "markerEvery": 1, 'markersize':5, 'label': 'FEM-22'}),
    # (fem_23,{"method": "FEM", "F":"Fy", "t":"t", "dt": 0.0015, "lineStyle": '--', "markerEvery": 1, 'markersize':5, 'label': 'FEM-23'}),
    # (fem_25,{"method": "FEM", "F":"Fy", "t":"t", "dt": 0.0015, "lineStyle": '--', "markerEvery": 1, 'markersize':5, 'label': 'FEM-25'}),
    # (fem_26,{"method": "FEM", "F":"Fy", "t":"t", "dt": 0.0015, "lineStyle": '--', "markerEvery": 1, 'markersize':5, 'label': 'FEM-26'}),
    # (fem_27,{"method": "FEM", "F":"Fy", "t":"t", "dt": 0.0015, "lineStyle": '--', "markerEvery": 1, 'markersize':5, 'label': 'FEM-27'}),


    # (fem_28,{"method": "FEM", "F":"Fy", "t":"t", "dt": 0.0015, "lineStyle": '--', "markerEvery": 1, 'markersize':5, 'label': 'FEM-28'}),
    # (fem_29,{"method": "FEM", "F":"Fy", "t":"t", "dt": 0.0015, "lineStyle": '--', "markerEvery": 1, 'markersize':5, 'label': 'FEM-29'}),
    # (fem_30,{"method": "FEM", "F":"Fy", "t":"t", "dt": 0.0015, "lineStyle": '--', "markerEvery": 1, 'markersize':5, 'label': 'FEM-30'}),
    # (fem_31,{"method": "FEM", "F":"Fy", "t":"t", "dt": 0.0015, "lineStyle": '--', "markerEvery": 1, 'markersize':5, 'label': 'FEM-31'}),
    # (fem_32,{"method": "FEM", "F":"Fy", "t":"t", "dt": 0.0015, "lineStyle": '--', "markerEvery": 1, 'markersize':5, 'label': 'FEM-32'}),
    # (fem_33,{"method": "FEM", "F":"Fy", "t":"t", "dt": 0.0015, "lineStyle": '--', "markerEvery": 1, 'markersize':5, 'label': 'FEM-33'}),
    # (fem_35,{"method": "FEM", "F":"Fy", "t":"t", "dt": 0.0015, "lineStyle": '--', "markerEvery": 1, 'markersize':5, 'label': 'FEM-35'}),
    # (fem_36,{"method": "FEM", "F":"Fy", "t":"t", "dt": 0.0015, "lineStyle": '--', "markerEvery": 1, 'markersize':5, 'label': 'FEM-36'}),
    # (fem_37,{"method": "FEM", "F":"Fy", "t":"t", "dt": 0.0015, "lineStyle": '--', "markerEvery": 1, 'markersize':5, 'label': 'FEM-37'}),
    # (fem_38,{"method": "FEM", "F":"Fy", "t":"t", "dt": 0.0015, "lineStyle": '--', "markerEvery": 1, 'markersize':5, 'label': 'FEM-38'}),
    # (fem_39,{"method": "FEM", "F":"Fy", "t":"t", "dt": 0.0015, "lineStyle": '--', "markerEvery": 1, 'markersize':5, 'label': 'FEM-39'}),
    # (fem_40,{"method": "FEM", "F":"Fy", "t":"t", "dt": 0.0015, "lineStyle": '--', "markerEvery": 1, 'markersize':5, 'label': 'FEM-40'}),




    # (isph_0,{"method": "SPH", "F":"Fz","t":"Time", "dt": 0.001,"lineStyle": '-', "markerEvery": 1, 'markersize':5, 'label': 'SPH-0'}),
    # (isph_1,{"method": "SPH", "F":"Fz","t":"Time", "dt": 0.001,"lineStyle": '-', "markerEvery": 1, 'markersize':5, 'label': 'SPH-1'}),
    #### (isph_2,{"method": "SPH", "F":"Fz","t":"Time", "dt": 0.001,"lineStyle": '-', "markerEvery": 1, 'markersize':5, 'label': 'SPH-2'}),
    # (isph_3,{"method": "SPH", "F":"Fz","t":"Time", "dt": 0.001,"lineStyle": '-', "markerEvery": 1, 'markersize':5, 'label': 'SPH-3'}),
    #### (isph_4,{"method": "SPH", "F":"Fz","t":"Time", "dt": 0.001,"lineStyle": '-', "markerEvery": 1, 'markersize':5, 'label': 'SPH-4'}),
    # (isph_5,{"method": "SPH", "F":"Fz","t":"Time", "dt": 0.001,"lineStyle": '-', "markerEvery": 1, 'markersize':5, 'label': 'SPH-5'}),
    #### (isph_6,{"method": "SPH", "F":"Fz","t":"Time", "dt": 0.001,"lineStyle": '-', "markerEvery": 1, 'markersize':5, 'label': 'SPH-6'}),
    #### (isph_7,{"method": "SPH", "F":"Fz","t":"Time", "dt": 0.001,"lineStyle": '-', "markerEvery": 1, 'markersize':5, 'label': 'SPH-7'}),
    ####(isph_8,{"method": "SPH", "F":"Fz","t":"Time", "dt": 0.001,"lineStyle": '-', "markerEvery": 1, 'markersize':5, 'label': 'SPH-8'}),
    # (isph_9,{"method": "SPH", "F":"Fz","t":"Time", "dt": 0.001,"lineStyle": '-', "markerEvery": 1, 'markersize':5, 'label': 'SPH-9'}),
    #### isph_10,{"method": "SPH", "F":"Fz","t":"Time", "dt": 0.001,"lineStyle": '-', "markerEvery": 1, 'markersize':5, 'label': 'SPH-10'}),
    #### isph_11,{"method": "SPH", "F":"Fz","t":"Time", "dt": 0.001,"lineStyle": '-', "markerEvery": 1, 'markersize':5, 'label': 'SPH-11'}),
    # (isph_12,{"method": "SPH", "F":"Fz","t":"Time", "dt": 0.001,"lineStyle": '-', "markerEvery": 1, 'markersize':5, 'label': 'SPH-12'}),
    # (isph_13,{"method": "SPH", "F":"Fz","t":"Time", "dt": 0.001,"lineStyle": '-', "markerEvery": 1, 'markersize':5, 'label': 'SPH-13'}),
    ####(isph_14,{"method": "SPH", "F":"Fz","t":"Time", "dt": 0.001,"lineStyle": '-', "markerEvery": 1, 'markersize':5, 'label': 'SPH-14'}),
    #### (isph_15,{"method": "SPH", "F":"Fz","t":"Time", "dt": 0.001,"lineStyle": '-', "markerEvery": 1, 'markersize':5, 'label': 'SPH-15'}),
    # (isph_16,{"method": "SPH", "F":"Fz","t":"Time", "dt": 0.001,"lineStyle": '-', "markerEvery": 1, 'markersize':5, 'label': 'SPH-16'}),
    ####(isph_17,{"method": "SPH", "F":"Fz","t":"Time", "dt": 0.001,"lineStyle": '-', "markerEvery": 1, 'markersize':5, 'label': 'SPH-17'}),
    ####(isph_18,{"method": "SPH", "F":"Fz","t":"Time", "dt": 0.001,"lineStyle": '-', "markerEvery": 1, 'markersize':5, 'label': 'SPH-18'}),
    #### (isph_19,{"method": "SPH", "F":"Fz","t":"Time", "dt": 0.001,"lineStyle": '-', "markerEvery": 1, 'markersize':5, 'label': 'SPH-19'}),
    # (isph_20,{"method": "SPH", "F":"Fz","t":"Time", "dt": 0.001,"lineStyle": '-', "markerEvery": 1, 'markersize':5, 'label': 'SPH-20'}),
    # (isph_21,{"method": "SPH", "F":"Fz","t":"Time", "dt": 0.001,"lineStyle": '-', "markerEvery": 1, 'markersize':5, 'label': 'SPH-21'}),
    # (isph_22,{"method": "SPH", "F":"Fz","t":"Time", "dt": 0.001,"lineStyle": '-', "markerEvery": 1, 'markersize':5, 'label': 'SPH-22'}),
    # (isph_23,{"method": "SPH", "F":"Fz","t":"Time", "dt": 0.001,"lineStyle": '-', "markerEvery": 1, 'markersize':5, 'label': 'SPH-23'}),
    ####( (isph_24,{"method": "SPH", "F":"Fz","t":"Time", "dt": 0.001,"lineStyle": '-', "markerEvery": 1, 'markersize':5, 'label': 'SPH-24'}),
    ####((isph_25,{"method": "SPH", "F":"Fz","t":"Time", "dt": 0.001,"lineStyle": '-', "markerEvery": 1, 'markersize':5, 'label': 'SPH-25'}),
    # (isph_26,{"method": "SPH", "F":"Fz","t":"Time", "dt": 0.001,"lineStyle": '-', "markerEvery": 1, 'markersize':5, 'label': 'SPH-26'}),
    #(isph_27,{"method": "SPH", "F":"Fz","t":"Time", "dt": 0.001,"lineStyle": '-', "markerEvery": 1, 'markersize':5, 'label': 'SPH-27'}),
    ####((isph_28,{"method": "SPH", "F":"Fz","t":"Time", "dt": 0.001,"lineStyle": '-', "markerEvery": 1, 'markersize':5, 'label': 'SPH-28'}),

])

fig = plt.figure(num=None, figsize=(16, 12),
                facecolor='w', edgecolor='k')
ax = fig.add_subplot(111)
# ax.set_title("")
ax.set_ylabel(r'$F_y (N)$')
ax.set_xlabel(r't($s$)')
ax.grid(color='k', linestyle='--', linewidth=0.2)
ax.autoscale(enable=True, axis='kx', tight=True)


for thisFile in files:
    data = pd.read_csv(thisFile)
    Force =data[files[thisFile]["F"]]
    if files[thisFile]['method']=="FEM":
        time=data[files[thisFile]["t"]]
        Force=FEMFunction(Force)
    if files[thisFile]['method']=="SPH":
        time=data[files[thisFile]["t"]][::2]
        Force=FEMFunction(Force[::2])


    plt.plot(time, Force,
              files[thisFile]["lineStyle"],
              linewidth=1.0, markersize=files[thisFile]["markersize"], markevery=files[thisFile]["markerEvery"],
              label=files[thisFile]['label'])



plt.legend( fancybox=True, shadow=True, ncol=1)

plt.savefig('Figure_Cylinder_FSI.png', bbox_inches='tight')
plt.show()
