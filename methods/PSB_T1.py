from dev_V2.six_qubit_QC_v2.system import six_dot_sample
from core_tools.sweeps.pulse_lib_wrappers.PSB_exp import run_qubit_exp
import pulse_lib.segments.utility.looping as lp
from core_tools.utility.variable_mgr.var_mgr import variable_mgr
from core_tools.sweeps.sweeps import scan_generic
from good_morning.calibrations.ultility import readout_convertor
from good_morning.fittings.fit_resonance import fit_resonance
from pulse_templates.coherent_control.single_qubit_gates.single_qubit_gates import single_qubit_gate_spec
from core_tools.data.ds.data_set import load_by_uuid

from dev_V2.six_qubit_QC_v2.VAR import variables

import numpy as np
import qcodes as qc

def PSB12_T1(state):
    '''
    Parameters
    ----------

    Returns
    -------
    function
        function to run the experiment

    '''
    s = six_dot_sample(qc.Station.default.pulse)
    
    gates, _311113, ST_anti_12, ST_anti_12_tc_high, ST_anti_56, ST_anti_56_tc_high, vSD1_threshold, vSD2_threshold = variables()
    anticrossing = list(ST_anti_12)
    # anticrossing = list(ST_anti_12)
    gates = list(gates)
    
    
    s.add(s.init12_FPGA)
    s.sequencer.add_ramp(gates, 100, _311113, ST_anti_12)
    
    if state is 'T_dd':
        s.add(s.q1.X180)
    
        
    if state is 'T_uu':
        s.add(s.q2.X180)
    
    
    waiting_time = lp.linspace(2e6,2e2,25, name='time', unit ='ns',axis=0)
    s.add_long_block(gates, waiting_time, anticrossing, sample_rate=1e9)

    

    
    loops = lp.linspace(1,15,15, name='looping', unit ='#',axis=1)
    s.add(s.wait(100+loops))
    s.add(s.read12)
    s.n_rep = 200
    
    return run_qubit_exp(f'PSB12 {state} T1 decay', s.sequencer)


def PSB56_T1(state):
    '''
    Parameters
    ----------

    Returns
    -------
    function
        function to run the experiment

    '''
    
    s = six_dot_sample(qc.Station.default.pulse)
    
    gates, _311113, ST_anti_12, ST_anti_12_tc_high, ST_anti_56, ST_anti_56_tc_high, vSD1_threshold, vSD2_threshold = variables()
    anticrossing = list(ST_anti_56)
    gates = list(gates)
    
    
    s.add(s.init56_FPGA)
    s.sequencer.add_ramp(gates, 100, _311113, ST_anti_56)
    
    if state is 'T_dd':
        s.add(s.q6.X180)
        
    if state is 'T_uu':
        s.add(s.q5.X180)
    
    waiting_time = lp.linspace(2e6,2e3,25, name='time', unit ='ns',axis=0)
    s.add_long_block(gates, waiting_time, anticrossing, sample_rate=1e9)

    

    
    loops = lp.linspace(1,20,20, name='looping', unit ='#',axis=1)
    s.add(s.wait(100+loops))
    s.add(s.read56)
    s.n_rep = 200
    
    return run_qubit_exp(f'PSB56 {state} T1 decay', s.sequencer)

#%%
from lmfit import Parameters, minimize
from lmfit.printfuncs import report_fit
import matplotlib.pyplot as plt
def Decay_formula(pars, x, data=None):
    amp = pars['amp'].value
    off = pars['off'].value
    T1 = pars['T1'].value
    
    model = amp*np.exp(-x/T1) + off

    if data is None:
        return model

    return np.nan_to_num((model - data)*1e6)
    
def fit_decay_T1(x,y, plot=False):
    fit_params = Parameters()
    fit_params.add('amp', value=0.4, max=1.0, min=-1.0)
    fit_params.add('off', value=0.1, max=1.0, min=-1.0)
    fit_params.add('T1', value=2e-5)
    
    out = minimize(Decay_formula, fit_params, args=(x,), kws={'data': y})
    if plot == True:
        fit = Decay_formula(out.params, x)    
        plt.plot(x,y,'o')
        plt.plot(x,fit,linewidth=5, alpha=0.5)
        plt.show()
        
    print(f'T1 = {out.params["T1"].value} s')
# %% Define and run experiment

if __name__ == "__main__":
    sequence, minstr, name = PSB12_T1('T_uu')
    ds = scan_generic(sequence, minstr, name=name).run()
    
    # sequence, minstr, name = PSB56_T1('T_dd')
    # ds = scan_generic(sequence, minstr, name=name).run()
    
    
    # ds_avg = ds['read12'].average('x')
    # fit_decay_T1(1e-9*ds_avg.x()[5:],ds_avg.y()[5:],True)
