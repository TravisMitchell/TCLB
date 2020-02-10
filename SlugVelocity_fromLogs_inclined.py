import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
def movingaverage(interval, window_size):
    window= np.ones(int(window_size))/float(window_size)
    return np.convolve(interval, window, 'same')
mavg = 120

if len(sys.argv) >= 2:
    myFile = sys.argv[1:]
else: 
    myFile = [r'Eo100\AnnularTaylorBubble_prmSweep_100_1_128_12000_65_Log_P00_00000010.csv']   # sys.argv[1]

endhit = 1

#endhit = 4000 # Eo10_Mo1e-4
#endhit = 5000 # Eo10_Mo1e-6
#endhit = 5600 # Eo20_Mo1e-6
#endhit = 5000 # Eo20_Mo1e-4
#endhit = 2500 # Eo20_Mo1e-2
#endhit = 7400 # Eo100_Mo1e-6
#endhit = 7800 # Eo400_Mo1e-4
#endhit = 7000 # Eo400_Mo1e-6
#endhit = 7200 # Eo700_Mo0.01
#endhit = 5500 # Eo200_Mo1e-6_15deg
print(myFile, len(myFile))

for i in myFile:
    data = pd.read_csv( i )
    
    its = data['Iteration'].values[:-endhit]
    try:
        bubFront = data['HeightCut'].values[:-endhit]
        print([len(bubFront), its[-1]])
    except:
        bubFront = data['BubbleFront'].values[:-endhit]
        print([len(bubFront), its[-1]])
    dPosition = bubFront[1:] - bubFront[:-1]
    instant_V = dPosition/(its[2]-its[1])
    instant_V_av = movingaverage( instant_V, mavg )

    average_V = data['GasTotalVelocityY'].values[:-endhit] / data['GasCells'].values[:-endhit]
    average_V_av = movingaverage( average_V, mavg )
    
    axs1 = plt.subplot(3,1,1)
    axs1.set_title(i[10:-4])
    axs1.plot(its[1:], instant_V);     axs1.plot(its[1:-1*int(mavg/2)], instant_V_av[:-1*int(mavg/2)],'r-')
    axs1.set_ylabel('velocity [lu/lt]');    
    axs1.text(max(its)/2, 0.25*max(instant_V_av), 'vel = ' + str(round(np.mean(instant_V_av[-1*int(mavg):-1*int(mavg/2)]),8)))
    axs2 = plt.subplot(3,1,2)
    axs2.plot(its, average_V);     axs2.plot(its[:-1*int(mavg/2)], average_V_av[:-1*int(mavg/2)],'r-')
    axs2.text(max(its)/2, 0.25*max(average_V_av), 'vel = ' + str(round(np.mean(average_V_av[-1*int(mavg):-1*int(mavg/2)]),8)))
    axs2.set_ylabel('velocity [lu/lt]');    
    axs3 = plt.subplot(3,1,3)
    axs3.plot(its[1:-1*int(mavg/2)], instant_V_av[:-1*int(mavg/2)],'b-')
    axs3.plot(its[:-1*int(mavg/2)], average_V_av[:-1*int(mavg/2)],'r-')
    axs3.set_ylabel('velocity [lu/lt]'); 
    axs3.set_xlabel('Iterations')
    plt.savefig(i[:-4] + '.png',format='png')
    plt.show()




