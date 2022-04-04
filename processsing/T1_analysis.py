import sys
sys.path.insert(1,r'c:\V2_code\dev_V2')
#%%
from core_tools.data.ds.data_set import load_by_uuid
from heating_uri.methods.T1_charac import  Decay_formula, fit_decay_T1
import matplotlib.pyplot as plt
import numpy as np
plot=False
temps = [15,100,200,250,300,350,400,500,600]
ids = [[1648584095408974076,1648587930020974076],
       [1648570439197974076,1648574270342974076],
       [1648663308962974076,1648667143390974076],
       [1648848295653974076,1648850450149974076],
       [1648809365891974076,1648813202027974076],
       [1648896673144974076,1648898472612974076],
       [1648905931271974076,1648907729851974076],
       [1648915747692974076,1648916327537974076],
       [1648917703940974076,1648920814892974076]]


dic_data = {'q1':np.empty((len(temps),2)),'q2':np.empty((len(temps),2))}
for ii, temp  in enumerate(temps):
    dss = ids[ii]
    for qubit, jj in enumerate(dss):
        qubit +=1
        ds = load_by_uuid(jj)
        ds_avg = ds[f'read{qubit}'].average('x')
        #it is 1-y for T1 and only y for T1_PSB
        params = fit_decay_T1(1e-9*ds_avg.x(), 1-ds_avg.y(), jj,plot=plot)
        dic_data[f'q{qubit}'][ii][0] = params['T1'].value
        dic_data[f'q{qubit}'][ii][1] = params['T1'].stderr
        

plt.figure()
for qubit in list(dic_data.keys()):
    plt.errorbar(temps,np.array(dic_data[qubit][:,0]), yerr=np.array(dic_data[qubit][:,1]),fmt='o-',label=qubit)
    plt.xlabel('T_mc [mK]')
    plt.ylabel('T1 [s]')
    
    # x = np.linspace(300,600,50)
    # plt.loglog(x, 5e14*x**-5,alpha= 0.5, linewidth= 5)
    # x = np.linspace(15,300,50)
    # plt.loglog(x, 7e2*x**-0.4, alpha= 0.5, linewidth= 5)
plt.legend()
plt.show()        