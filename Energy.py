# Packages
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import gridspec

# Energy data extraction
data=np.loadtxt('energy.txt')

t=data[:,0]
u=data[:,1]
P=data[:,2]
K=data[:,3]
Etot=u+P+K
DE=(Etot-Etot[0])/abs(Etot[0])*100



fig = plt.figure()
# Set height ratios for subplots
gs = gridspec.GridSpec(2, 1, height_ratios=[2, 1]) 

# The first subplot
ax0 = plt.subplot(gs[0])

line0, = ax0.plot(t, u, color='b',linestyle='-.')
line1, = ax0.plot(t, P, color='r',linestyle='-.')
line2, = ax0.plot(t, K, color='g',linestyle='-.')
line3, = ax0.plot(t, Etot, color='k',linestyle='-.')
plt.xlabel('Time')
plt.ylabel('Energy')
plt.grid()


# The second subplot
# Shared axis X
ax1 = plt.subplot(gs[1], sharex = ax0)
line4, = ax1.plot(t, DE, color='k', linestyle='-.')
plt.setp(ax0.get_xticklabels(), visible=False)
# Remove last tick label for the second subplot
yticks = ax1.yaxis.get_major_ticks()
yticks[-1].label1.set_visible(False)

# Put legend on first subplot
ax0.legend((line0, line1, line2, line3), ('Internal', 'Potential','Kinetic','Total'), loc='lower left')
plt.xlabel('Time')
plt.ylabel(r'$\Delta E_{tot}/|E_{tot,0}|$x100')
plt.xlim(0,3)
# Remove vertical gap between subplots
plt.subplots_adjust(hspace=.0)
plt.grid()
#plt.savefig('Energys')
plt.show()
