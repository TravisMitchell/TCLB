import sys
import h5py
import numpy as np
import math
import matplotlib.pyplot as plt
import pandas as pd

if len(sys.argv) >= 2:
    myFile = sys.argv[1:]
else: 
    myFile = r'AnnularTaylorBubble_prmSweep_10_1e-06_128_12000_HDF5_01200010.h5'   # sys.argv[1]

for i in range(len(myFile)):
    data = h5py.File(myFile[i],'r')
#    print(myFile[i])
    PF = np.array(data['PhaseField'])
    U  = np.array(data['U'])
    Wall = np.array(data['WallBoundary'])

#    print('Starting integration...')
    integral_U  = np.sqrt(    np.sum(np.sum(U[:,:,:,0], 2), 0 ) **2 \
			    + np.sum(np.sum(U[:,:,:,1], 2), 0 ) **2 \
			    + np.sum(np.sum(U[:,:,:,2], 2), 0 ) **2 )
    integral_pf = np.sum(np.sum(PF[:,:,:], 2), 0 )
#    print('Finished Integration import...')
    
    if max(integral_U) > 0 and max(integral_pf) > 0:
        integral_U = integral_U/max(integral_U)
        integral_pf= integral_pf/max(integral_pf)
    else:
        continue

    # Take search of interface at 0.95
    myLim = 0.95
    n = len(integral_pf)
    upper_Bub = n
    lower_Bub = n
    s_u = 0.0

    for j in range(1,n):
        if integral_pf[n-j-1] < myLim and integral_pf[n-j] > myLim and s_u == 0.0:
            upper_Bub = n-j-1
            s_u = 1.0
        if integral_pf[n-j-1] > myLim and integral_pf[n-j] < myLim:
            lower_Bub = n-j-1
            break
    # Take wake at 0.05 
    myLim_U = 0.05
    wake_end = lower_Bub
    for j in range(lower_Bub):
        if integral_U[lower_Bub - j] < myLim_U:
            wake_end = lower_Bub - j
            break

#    print("Bubble front %f and bubble end %f so bubble length %f" % (upper_Bub, lower_Bub, upper_Bub - lower_Bub))
#    print("Wake length %f, characteristic length %f so dimensionless wake = %f" % (wake_end, U.shape[0], wake_end/float(U.shape[0])))
    print("%s, bubble front %.1f, len(bub) %.1f, wake length %.1f=%f" % (myFile[i][38:], upper_Bub, upper_Bub-lower_Bub, (lower_Bub-wake_end), (lower_Bub-wake_end)/128.0) )
    plt.figure(1)
    plt.clf()
    plt.plot(range(len(integral_U)), integral_U, 'r')
    plt.plot(range(len(integral_pf)), integral_pf, 'b')
    plt.plot([lower_Bub, upper_Bub, wake_end], [myLim, myLim, myLim_U], 'k.')
    plt.savefig(myFile[i] + '_int.png')

    df = pd.DataFrame({'PF_int': integral_pf, 'U_int':integral_U})
    out = myFile[i] + '.csv'
    df.to_csv(out)

