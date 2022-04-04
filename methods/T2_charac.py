from dev_V2.six_qubit_QC_v2.system import six_dot_sample
from core_tools.sweeps.pulse_lib_wrappers.PSB_exp import run_qubit_exp
import pulse_lib.segments.utility.looping as lp



from core_tools.sweeps.sweeps import scan_generic, do0D, do1D
from pulse_templates.coherent_control.single_qubit_gates.single_qubit_gates import single_qubit_gate_spec
from pulse_templates.oper.operators import jump, wait


import qcodes as qc
import numpy as np


def T2_ramsey(seg, gate_set, t_wait, f_bg_oscillations):
    '''
    perform T2* measurement on the qubit

    Args:
        seg (segment_container) : segment container
        gate_set (single_qubit_gate_spec) : object containg instruction of the gate set
        t_wait (double) : time to wait
        f_bg_oscillations (double) : freq at which the the qubit needs to oscilate 
    '''
    gate_set.X.build(seg)
    gate_set.wait(seg, t_wait)
    gate_set.Z(t_wait*1e-9*f_bg_oscillations*np.pi*2).build(seg)
    gate_set.X.build(seg)
    

def T2star(target_q):
    s = six_dot_sample(qc.Station.default.pulse)
    s.sequencer.n_rep = 1000
    
    s.init(target_q)
   
    s.add_func(T2_ramsey, s[f'q{target_q}'], lp.linspace(6e3,10,75, name='time', unit ='ns',axis=0),1e6)  
    
    loops = lp.linspace(1,20,20, name='looping', unit ='#',axis=1)
    s.add(s.wait(100e3+loops))
    s.read(target_q)
    
    
    return run_qubit_exp(f'Ramsey Qubit {target_q}', s.sequencer)




def T2_CPMG_t_tot(seg, gate_set, t_wait, N_rep,f_bg_oscillations):
    '''
    return a CPMG sequence given tot total waiting time

    Args:
        seg (segment_container) : segment container
        gate_set (single_qubit_gate_spec) : object containg instruction of the gate set
        t_wait (double) : total time waited in the pulse sequence
        N_rep (double) : amount of X gates you want to do
    '''
    
    final_wait_time =  (t_wait-2*gate_set.X.t_pulse)/N_rep/2 - gate_set.X.t_pulse
    # final_wait_time  = t_wait
    gate_set.X.build(seg)
    for i in range(N_rep):
        gate_set.wait(seg,final_wait_time)
        gate_set.Y2.build(seg)
        gate_set.wait(seg, final_wait_time)
    # gate_set.Z(final_wait_time*1e-9*f_bg_oscillations*np.pi*2).build(seg)
    gate_set.X.build(seg)


def T2CPMGn(target_q, N_rep,ini_t,end_t):
    s = six_dot_sample(qc.Station.default.pulse)

    s.init(target_q)
    
    s.add_func(T2_CPMG_t_tot, s[f'q{target_q}'], lp.geomspace(end_t,ini_t,60, name='time', unit ='ns',axis=0), N_rep, 1e6)

    loops = lp.linspace(1,20,20, name='looping', unit ='#',axis=1)
    s.add(s.wait(100e3+loops))

    s.read(target_q)
    s.sequencer.n_rep = 500
    
    return run_qubit_exp(f'CPMG Qubit {target_q}, N_rep={N_rep}', s.sequencer)

# %%

    
#%%
from lmfit import Parameters, minimize
from lmfit.printfuncs import report_fit
import matplotlib.pyplot as plt
def Ramsey_formula(pars, x, data=None):
    amp = pars['amp'].value
    off = pars['off'].value
    omega = pars['omega'].value
    phi = pars['phi'].value
    T2 = pars['T2'].value

    model = amp*np.cos(2*np.pi*omega*x + phi)*np.exp(-(x/T2)**2) + off

    if data is None:
        return model

    return np.nan_to_num((model - data)*1e6)
    
def fit_Ramsey(x,y, title='', plot=False):
    fit_params = Parameters()
    fit_params.add('amp', value=0.4, max=1.0, min=-1.0)
    fit_params.add('off', value=0.4, max=1.0, min=-10.0)
    fit_params.add('omega', value=2.2e5)
    fit_params.add('phi', value=0.0)
    fit_params.add('T2', value=2e-6, min= 0)
    
    out = minimize(Ramsey_formula, fit_params, args=(x,), kws={'data': y})
    if plot == True:
        plt.figure()
        fit = Ramsey_formula(out.params, x)    
        plt.plot(x,y,'o')
        plt.plot(x,fit,linewidth=5, alpha=0.5)
        plt.title(title)
        plt.show()
    print(f'\nT2* = {out.params["T2"].value}')
    print(f'Omega = {out.params["omega"].value}')
    return out.params
    
    
def Hahn_decay(pars, x, data=None):
    amp = pars['amp'].value
    off = pars['off'].value
    T2 = pars['T2'].value
    n = pars['n'].value

    model = off + amp*np.exp(-(x/T2)**n)
    if data is None:
        return model

    return (model - data)*1e6
    
def fit_Hahn(x,y, plot=False):
    fit_params = Parameters()
    fit_params.add('amp', value=0.3, max=1.0, min=-1.0)
    fit_params.add('off', value=0.5, max=1.0, min=-1.0)
    fit_params.add('T2', value=1e-4)
    fit_params.add('n', value=2)

    out = minimize(Hahn_decay, fit_params, args=(x,), kws={'data': y})
    if plot == True:
        fit = Hahn_decay(out.params, x)    
        plt.semilogx(x,y, 'o')
        plt.semilogx(x,fit,linewidth=5, alpha=0.5)
        plt.semilogx(x, off + amp*np.exp(-(x/T2)**n),linewidth=5, alpha=0.5)
        plt.show()
    # print(out.params)
    print(f'T2 = {out.params["T2"].value*1e6} us')

#%%
if __name__ == "__main__":
    from core_tools.sweeps.sweeps import scan_generic, do0D, do1D
    #Ramsey
    i=2
    # sequence, minstr, name = T2star(i)
    # ds = scan_generic(sequence, minstr, name=name).run() 
    ds_avg = ds[f'read{i}'].average('x')
    fit_Ramsey(ds_avg.x()*1e-9, ds_avg.y(), True)
    
        
    
    #CPMG
    N_rep = np.power(2,np.array(range(5,3)), dtype=int)
    for jj,reps in enumerate(N_rep):
        ini_t=10e2+reps*2*10e2
        # end_t = 10e5+reps*2.5*10e5
        # ini_t =10e2
        end_t = 9e7
        sequence, minstr, name = T2CPMGn(i, reps,ini_t,end_t)
        ds = scan_generic(sequence, minstr, name=name).run()
        ds_avg = ds[f'read{i}'].average('x')
        fit_Hahn(ds_avg.x()*1e-9, ds_avg.y(), True)
        
    

