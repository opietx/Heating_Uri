import sys
sys.path.insert(1,r'c:\V2_code\dev_V2')
#%%
from core_tools.data.ds.data_set import load_by_uuid
from heating_uri.methods.T2_charac import  Ramsey_formula, fit_Ramsey
import matplotlib.pyplot as plt
import numpy as np
temps = [15,100,200,250,300]
ids = [[1647549612187974076,1647550547152974076],
       [1648557419650974076,1648557900168974076],
       [1648484854195974076,1648485316244974076],
       [1648852264799974076,1648852735723974076],
       [1648821042156974076,1648822592211974076]]


dic_data = {'q1':np.empty(len(temps)),'q2':np.empty(len(temps))}
for ii, temp  in enumerate(temps):
    dss = ids[ii]
    for qubit, jj in enumerate(dss):
        qubit +=1
        print(dss,qubit)
        ds = load_by_uuid(jj)
        ds_avg = ds[f'read{qubit}'].average('x')
        params = fit_Ramsey(ds_avg.x()*1e-9, ds_avg.y(),jj, True)
        dic_data[f'q{qubit}'][ii] = params['T2'].value
        

plt.figure()
for qubit in list(dic_data.keys()):
    plt.plot(temps,np.array(dic_data[qubit])*1e6,'o',label=qubit)
    plt.xlabel('T_mc [mK]')
    plt.ylabel('T2* [us]')
plt.legend()
plt.show()     