from dev_V2.six_qubit_QC_v2.system import six_dot_sample
from core_tools.sweeps.pulse_lib_wrappers.PSB_exp import run_qubit_exp
import pulse_lib.segments.utility.looping as lp
from pulse_templates.coherent_control.RB_single.RB_seg_generator import generature_single_qubit_RB

import qcodes as qc
import numpy as np
def RB_run(target_q):
    s = six_dot_sample(qc.Station.default.pulse)

    s.add(s.init12_FPGA)
    

    s.add_func(generature_single_qubit_RB, s[f'q{target_q}'], lp.linspace(2,300,100, name='nth set', unit ='#'), 10, 'XY')
    s.add(s.wait(400e3))

    # s.add(s.read123_FPGA_high_fid)
    s.read(target_q)
    return run_qubit_exp(f'RB Qubit {target_q}', s.sequencer)



#%%
from lmfit import Parameters, minimize
from lmfit.printfuncs import report_fit
import matplotlib.pyplot as plt
def rb_formula(pars,x, data=None):
    A = pars['A']
    B = pars['B']
    F = pars['Fid']
    
    model = A*(2*F-1)**x + B
    if data is None:
        return model
    
    return model-data
    
def fit_RB(x,y, plot=False):
    fit_params = Parameters()

    fit_params.add('A', value=-0.5, max=0.5) # read out fidelity/2
    fit_params.add('B', value=0.5, min=0.4,max=0.6) # offset of the visibility
    fit_params.add('Fid', value=0.99)

    out = minimize(rb_formula, fit_params, args=(x,), kws={'data': y})
    
    if plot == True:
        fit = rb_formula(out.params, x)
        plt.semilogx(x,y,'o')
        plt.semilogx(x,fit)
        plt.show()
    print(out.params)
    print(f'fidelity = {out.params["Fid"].value*100}%')
    print(f'single gate fidelity = {(1-(1-out.params["Fid"].value)/1.875)*100}%')

#%% RUn
if __name__ == "__main__":    
    from core_tools.sweeps.sweeps import scan_generic, do0D, do1D
    
    for ii in [1,2]:
        sequence, minstr, name = RB_run(ii)
        ds = scan_generic(sequence, minstr, name=name).put()
        
    
    ds_avg = ds['read1'].average('x')
    
    fit_RB(ds_avg.x()[:-2], ds_avg.y()[:-2], True)
