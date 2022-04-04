from dev_V2.six_qubit_QC_v2.system import six_dot_sample
from core_tools.sweeps.pulse_lib_wrappers.PSB_exp import run_qubit_exp
import pulse_lib.segments.utility.looping as lp
from pulse_templates.coherent_control.single_qubit_gates.single_qubit_gates import single_qubit_gate_spec




import qcodes as qc
import numpy as np


def T1(target_q, t_final=1e9):
    s = six_dot_sample(qc.Station.default.pulse)
    s.n_rep = 500
   
    
        
    if target_q in [1,2]: 
        s.add(s.init12_FPGA)
    else:
        s.add(s.init56_FPGA)
    
    if  target_q == 2: 
        s.q1.X180
        s.q2.X180
    if  target_q == 5: 
        s.q5.X180
        s.q6.X180

    # s.add(s.q5.X180)   
    waiting_time = lp.linspace(21e3,t_final,15, name='time', unit ='ns',axis=0)
    s.wait_long(waiting_time, sample_rate=1e8)
    
    # waiting_time = lp.linspace(210e3,t_final,10, name='time', unit ='ns',axis=0)
    # s.wait_long(waiting_time, sample_rate=1e7)
    
    
    loops = lp.linspace(1,10,10, name='looping', unit ='#',axis=1)
    s.add(s.wait(10+loops))
    
    s.read(target_q)   

    s.sequencer.reset_time()
    return run_qubit_exp(f'Decay T1 Qubit {target_q}', s.sequencer)


    
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
    
def fit_decay_T1(x,y, title='',plot=False):
    fit_params = Parameters()
    fit_params.add('amp', value=0.4, max=1.0, min=0.0)
    fit_params.add('off', value=0.7, max=1.0, min=0.0)
    fit_params.add('T1', value=x[int(len(x)/2)])
    
    out = minimize(Decay_formula, fit_params, args=(x,), kws={'data': y})
    if plot == True:
        fit = Decay_formula(out.params, x)    
        plt.figure()
        plt.plot(x,y,'o')
        plt.plot(x,fit,linewidth=5, alpha=0.5)
        plt.title(title)
        plt.show()
        
    print(f'T1 = {out.params["T1"].value} s')
    return out.params
    

#%%Run
if __name__ == "__main__":
    from core_tools.sweeps.sweeps import scan_generic, do0D, do1D
    for i in [1,2]:
        sequence, minstr, name = T1(i,800e6)
        ds = scan_generic(sequence, minstr, name=name).run()
 
    ds_avg = ds[f'read{i}'].average('x')
    a=fit_decay_T1(1e-9*ds_avg.x()[1:], 1-ds_avg.y()[1:], plot=True)
    
    
    nums = []
    for loop in range(len(ds['read1'].x())):
        ds_avg = ds['read1'][loop]
        a=fit_decay_T1(1e-9*ds_avg.x()[1:], 1-ds_avg.y()[1:], plot=True)
        nums.append([a['T1'].value,a['T1'].stderr])
        
        
    
    
    

