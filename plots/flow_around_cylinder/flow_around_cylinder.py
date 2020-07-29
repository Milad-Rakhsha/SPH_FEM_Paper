import scipy.signal as signal
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
matplotlib.rcParams.update({'font.size': 24})

def filterisph(x, N=5, Wn=0.99):
    B, A = signal.butter(N, Wn, output='ba')
    smooth = signal.filtfilt(B,A, x)
    return smooth

def filterfem(x,N=5, Wn=0.999):
    B, A = signal.butter(N, Wn, output='ba')
    smooth = signal.filtfilt(B,A, x)
    return smooth

isph_=("isph.csv.414176")
isph_=("isph.csv.409398")
isph_=("isph.csv.426921")
# isph_=("isph.csv.428327")
# isph_=("isph.csv.428287")
# isph_=("isph.csv.427031")
# isph_=("isph.csv.427011")
# isph_=("isph.csv.426991")
# isph_=("isph.csv.426985")
# isph_=("isph.csv.426973")
# isph_=("isph.csv.426969")
# isph_=("isph.csv.426958")
isph_=("isph.csv.426921")
# isph_=("isph.csv.426888")

# isph_=("isph.csv")
# isph_=("wcsph.csv")
fem_=("fem.csv")


fig = plt.figure(num=None, figsize=(16, 12),
                facecolor='w', edgecolor='k')
ax = fig.add_subplot(111)
ax.set_title("Flow Around Cylinder", fontsize=36)
ax.set_ylabel(r'$Cd=Fx/(1/2\rho U^2)$', fontsize=32)
ax.set_xlabel('time($s$)', fontsize=32)
ax.grid(color='k', linestyle='--', linewidth=0.5)
ax.autoscale(enable=True, axis='kx', tight=True)
ax.set(xlim=[0,10])

F_MAX=(0.5*100*0.95*0.95)
isph = pd.read_csv(isph_)
fem = pd.read_csv(fem_)/F_MAX

cd_isph=filterisph(isph["Fx"])/0.08/F_MAX
plt.plot(isph["Time"],cd_isph ,'-',
          linewidth=1.5, markersize=5, markevery=100,
          markerfacecolor='none', label='ISPH')

N=fem["Fx"].shape[0]
time=np.arange(N)*0.0005
plt.plot(time, filterfem(fem["Fx"]),'-o',
          linewidth=1.5, markersize=5, markevery=100,
          markerfacecolor='none', label='FEM')


plt.legend( fancybox=True, shadow=True, ncol=1)
plt.savefig('Figure_flow_around_cylinder.png', bbox_inches='tight')
plt.show()
