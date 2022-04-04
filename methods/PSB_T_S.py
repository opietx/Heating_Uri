from dev_V2.six_qubit_QC_v2.system import six_dot_sample
from core_tools.sweeps.pulse_lib_wrappers.PSB_exp import run_qubit_exp

import qcodes as qc
import numpy as np

def PSB12_S_T(state='S', readout_time = 10e3 , ramp_time=100, name=None):
    s = six_dot_sample(qc.Station.default.pulse)
    s.n_rep = 100000
    s.read12.measurement.t_measure = readout_time
    s.read12.measurement.ramp = ramp_time
    s.add(s.init12_FPGA)
    
    s.add(s.pre_pulse)
    if state is 'T_dd':
        s.add(s.q1.X180)
    if state is 'T_uu':
        s.add(s.q2.X180)
        
    # loops = lp.linspace(1,4,4, name='looping', unit ='#',axis=0)
    loops=0
    s.add(s.wait(10+loops))
    
    s.add(s.read12)
    return run_qubit_exp(f'PSB12 only {state}, thr= {name}', s.sequencer)


def PSB56_S_T(state='S', readout_time = 10e3 ,name=None):
    s = six_dot_sample(qc.Station.default.pulse)
    s.n_rep = 100000
    s.read56.measurement.t_measure = readout_time
    s.add(s.init56_FPGA)
    
    s.add(s.pre_pulse)
    if state is 'T_dd':
        s.add(s.q6.X180)
    if state is 'T_uu':
        s.add(s.q5.X180)
    
    s.add(s.wait(10e3))
    s.add(s.read56)
    return run_qubit_exp(f'PSB56 only {state}, thr= {name}', s.sequencer)



# %%


from core_tools.sweeps.sweeps import scan_generic, do0D, do1D
from core_tools.utility.variable_mgr.var_mgr import variable_mgr
import matplotlib.pyplot as plt
import numpy as np


def plot_hist(data):
    plt.figure()        
    for state in list(data.keys()):
        plt.hist(data.get(state), bins=int(np.sqrt(len(data.get(state)))), density=True, alpha=0.2, label=state)
    plt.legend()
        
    
    
if __name__ == "__main__":

    var_mngr = variable_mgr()
    thr1 = getattr(var_mngr, 'threshold_SD1')
    # thr2 = getattr(var_mngr, 'threshold_SD2')
    
    readout_time = 10e3
    
    
    states = ['S', 'T_uu', 'T_dd']
    data = {key: None for key in states}
    for state in states:
        # sequence, minstr, name = PSB56_S_T(state,readout_time=readout_time, name=round(thr2,2))
        sequence, minstr, name = PSB12_S_T(state,readout_time=readout_time, name=round(thr1,2))
        ds = scan_generic(sequence, minstr, name=name).run()




        pre_select_data =  ds['RAW init12_M2 (SD1_IQ:1)'].y()
        selection = np.where(pre_select_data< thr1) #this selects the proper initialization into Singlets 
        print('Selected {}'.format(round(len(selection[0])/len(pre_select_data)*100,2)))
        data_raw = ds['RAW read12 (SD1_IQ:2)'].y() # from the states that were regarded as good, we have these measurements
        data[state] = data_raw[selection]
    
    plot_hist(data)
    
    
    

