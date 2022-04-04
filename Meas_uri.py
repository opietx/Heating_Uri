import sys
sys.path.insert(1,r'c:\V2_code\dev_V2')

from pulse_templates.utility.plotting import plot_seg
#%% Load ds
from core_tools.data.ds.data_set import load_by_uuid
ds = load_by_uuid(1648663308962974076)


#%% T1
from heating_uri.methods.T1_charac import  *
from core_tools.sweeps.sweeps import scan_generic, do0D, do1D

qubits=[1,2]
for i in qubits:
    sequence, minstr, name = T1(i,105e6)
    ds = scan_generic(sequence, minstr, name=name).put()



#%% T2*
from heating_uri.methods.T2_charac import  *
from core_tools.sweeps.sweeps import scan_generic, do0D, do1D
plotting_enabled = True

qubits=[1,2]
for i in qubits:
    sequence, minstr, name = T2star(i)
    ds = scan_generic(sequence, minstr, name=name).run()



#%% CPMG
plotting_enabled = True
N_reps  = np.power(2,np.arange(4,-1,-1), dtype=int)

qubits=[1,2]
for i in qubits:
    for reps in N_reps:
        sequence, minstr, name = T2CPMGn(i, reps,6e5,6e3)
        ds = scan_generic(sequence, minstr, name=name).run()
    
    

#%% J symmetric point
from heating_uri.methods.J_char import  *

plotting_enabled = False
pairs = [12]
n_points = 10
start = 0.2
for pair in pairs:
    ds = J_symmmetry(pair, start=start, N_J = n_points, plot=plotting_enabled)


#%% PSB T1
from heating_uri.methods.PSB_T1 import  *
states = ['T_uu','T_dd']
for state in states:
    sequence, minstr, name = PSB12_T1(state)
    ds = scan_generic(sequence, minstr, name=name).run()

#for 56
# for state in states:
    # sequence, minstr, name = PSB56_T1(state)
    # ds = scan_generic(sequence, minstr, name=name).run()

#%%PSB analysis S-T-T
'''
Measure for 3 readout times,
for both sides
'''
from heating_uri.methods.PSB_T_S import  *
from core_tools.utility.variable_mgr.var_mgr import variable_mgr


var_mngr = variable_mgr()
states = ['S', 'T_uu', 'T_dd']
readout_time = 1e3


thr = getattr(var_mngr, 'threshold_SD1')    
for state in states:
    sequence, minstr, name = PSB12_S_T(state,readout_time=readout_time, name=round(thr,2))
    ds = scan_generic(sequence, minstr, name=name).run()


# thr = getattr(var_mngr, 'threshold_SD2')
# for state in states:
    # sequence, minstr, name = PSB56_S_T(state,readout_time, name=round(thr,2))
    # ds = scan_generic(sequence, minstr, name=name).run()




#%%Charge noise
'''If one want to measure charge noise must c'''
from heating_uri.methods.Charge_noise import  *


gate = station.gates.vSD2_P
v_start = gate()-7.5
v_stop =  gate()+7.5

# station.gates.V_src_1(100)


#To measure the DC coulom peaks
inst= [station.Idc_2]
ds = do1D(gate, v_start, v_stop, 200, 0.01, *inst,reset_param=True,name = f'1D gate scan of {gate.name}').run()

#To acquire charge noise
# measure_slope(gate, v_start, v_stop, station.keithley_3, 1)#lever arm of 1, it was 0.013 before

# station.gates.V_src_1(0)
