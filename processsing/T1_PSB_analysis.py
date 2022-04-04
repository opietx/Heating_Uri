import sys
sys.path.insert(1,r'c:\V2_code\dev_V2')
#%%
from core_tools.data.ds.data_set import load_by_uuid
from heating_uri.methods.T1_charac import  Decay_formula, fit_decay_T1
import matplotlib.pyplot as plt
import numpy as np
temps = [15,100,200,300]
ids = [[1648201848655974076,1648202272632974076],
       [1648563190746974076,1648563661289974076],
       [1648544773265974076,1648545706922974076],
       [1648828733548974076,1648829194873974076]]


dic_data = {'T_uu':np.empty(len(temps)),'T_dd':np.empty(len(temps))}
for ii, temp  in enumerate(temps):
    dss = ids[ii]
    for qubit, jj in enumerate(dss):
        ds = load_by_uuid(jj)
        ds_avg = ds['read12'].average('x')
        params = fit_decay_T1(1e-9*ds_avg.x(), ds_avg.y(),jj, True)
        dic_data[list(dic_data.keys())[qubit]][ii] = params['T1']
        

plt.figure()
for qubit in list(dic_data.keys()):
    plt.errorbar(temps,np.array(dic_data[qubit])*1e6,'o',label=qubit)
    plt.xlabel('T_mc [mK]')
    plt.ylabel('T1_PSB [us]')
plt.legend()
plt.show()

        