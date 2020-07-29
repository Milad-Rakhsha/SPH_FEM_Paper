from decimal import Decimal
import matplotlib
# matplotlib.use('TkAgg') 
matplotlib.use('Agg')  # Must be before importing matplotlib.pyplot or pylab!
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from collections import OrderedDict
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
matplotlib.rcParams.update({'font.size': 18})

    


isph_0=("isph.txt")
fem_0=("fem.txt")


files = OrderedDict([#
    (isph_0,{"method": "ISPH", "shift": +0.65, "dT": 0.05, "lineStyle": 'b-', "markerEvery": 1, 'markersize':5, 'label': 'ISPH'}),
    (fem_0,{"method": "FEM", "shift":-2.011,"dT": 0.02, "lineStyle": 'r-', "markerEvery": 1, 'markersize':5, 'label': 'FEM'}),

])

fig = plt.figure(num=None, figsize=(16, 12),
                facecolor='w', edgecolor='k')
ax = fig.add_subplot(111)
ax.set_title("Dam Break")
ax.set_ylabel(r'Front Position ($m$)')
ax.set_xlabel(r't ($s$)')
ax.grid(color='k', linestyle='--', linewidth=0.5)
ax.autoscale(enable=True, axis='kx', tight=True)


for thisFile in files:
    data = pd.read_csv(thisFile)
    # time=data["Time"]
    x =data["x"]+files[thisFile]["shift"]
    t=np.linspace(0,x.shape[0]+1,x.shape[0])*files[thisFile]["dT"]
    plt.plot(t, x,
              files[thisFile]["lineStyle"],
              linewidth=2, markersize=files[thisFile]["markersize"], markevery=files[thisFile]["markerEvery"],
              label=files[thisFile]['label'])
              


plt.legend( fancybox=True, shadow=True, ncol=1)
ax.set(xlim=[0, 0.8])

plt.savefig('Figure_damBreak.png', bbox_inches='tight')
plt.show()




