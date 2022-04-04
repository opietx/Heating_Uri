from good_morning.fittings.J_versus_voltage import fit_J, J_to_voltage, voltage_to_J
from good_morning.calibrations.ultility import get_target, readout_convertor
from good_morning.fittings.fit_rabi_osc import fit_ramsey

from core_tools.sweeps.pulse_lib_wrappers.PSB_exp import run_qubit_exp
from core_tools.utility.variable_mgr.var_mgr import variable_mgr
from core_tools.sweeps.sweeps import scan_generic, do0D, do1D
from core_tools.job_mgnt.job_mgmt import job_wrapper

from pulse_templates.coherent_control.two_qubit_gates.cphase import cphase_basic

from dev_V2.six_qubit_QC_v2.system import six_dot_sample
import pulse_lib.segments.utility.looping as lp

import qcodes as qc
import numpy as np

import good_morning.static.J12 as J12
import good_morning.static.J23 as J23
import good_morning.static.J34 as J34
import good_morning.static.J45 as J45
import good_morning.static.J56 as J56


@job_wrapper
def J_symmmetry(target, start=0.3,  N_J=5, plot=False):
	'''
	calibrates a single qubit phase for a target qubit

	Args:
		target (int) : qubit pair to target
	'''
	var_mgr = variable_mgr()
	s = six_dot_sample(qc.Station.default.pulse)

	target = str(int(target))
	s.init(target[0], target[1])

	s.add(s.wait(100))
	s.add(s.pre_pulse)
	s.add(s.wait(100))


	J_target = globals()[f'J{target}']
	barrier = lp.linspace(1, start, N_J, axis=1, name=f'vB{target[0]}', unit='mV')
	gate_time = lp.linspace(0, getattr(
	    var_mgr, f'cphase_time_{target}')*50, 400, axis=0, name='time', unit='ns')

	s.add(s[f'q{target[0]}'].X90)
	s.add_func(cphase_basic, J_target.gates, tuple(
	    [0]*len(J_target.gates)), J_target.barrier_perc_to_voltage(barrier), t_gate=gate_time/2, t_ramp=10)
	s.add(s[f'q{target[0]}'].X180)
	s.add(s[f'q{target[1]}'].X180)
	s.add_func(cphase_basic, J_target.gates, tuple(
	    [0]*len(J_target.gates)), J_target.barrier_perc_to_voltage(barrier), t_gate=gate_time/2, t_ramp=10)
	s.add(s[f'q{target[0]}'].X90)

	s.add(s.wait(100))

	s.read(target[0],target[1])

	sequence, minstr, name = run_qubit_exp(f'J on_off pair {target}', s.sequencer)
	ds = scan_generic(sequence, minstr, name=name).put()
# 	fit_J_evolution(target,ds,plot)
    
	return ds


	
def fit_J_evolution(target,ds, plot=False):
    
    time = 1e-9*ds(readout_convertor(f'read{target[0]}')).y()
    probabilities = ds(readout_convertor(f'read{target[0]}')).z()
    
	
    J_meas=[]
    for i in range(len(probabilities)):
        time_fit = time[5:]
        probabilities_fit = probabilities[i][5:]
        J_meas += [2/(fit_ramsey(time_fit, probabilities_fit, plot=plot)*2)*1e9]

    barrier_percentage = ds(readout_convertor(f'read{target[0]}')).x()
    mask  = np.where(np.array(J_meas) <15e6)
    J_V_off, J_max, alpha = fit_J(barrier_percentage[mask], np.array(J_meas)[mask], plot=True)
    


#%%

if __name__ == "__main__":
    pair = 12
    plotting_enabled = False
    ds = J_symmmetry(pair,start=0.1, N_J=9, plot=plotting_enabled)
    

   
