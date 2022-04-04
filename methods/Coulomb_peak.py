from core_tools.utility.variable_mgr.var_mgr import variable_mgr
from core_tools.job_mgnt.job_mgmt import job_wrapper

from core_tools.sweeps.sweeps import do0D

from core_tools.GUI.keysight_videomaps.data_getter.scan_generator_Keysight import  construct_1D_scan_fast

@job_wrapper
def SD1_peak(swing=20, verbose = False):
    var_mgr = variable_mgr()
    readout_point = getattr(var_mgr, 'PSB12_SD1')#get the calibrated readout point
    
    #add awg sweep #sampling rate??
    # fast_sweep = construct_1D_scan_fast('vSD1_P', 10, 100, 1e3, False, station.pulse, station.dig, [1], 1e8, 500e6, acquisition_delay_ns=0, enabled_markers=['M_SD1', 'M_SD2'])
    fast_sweep = construct_1D_scan_fast('vSD1_P', swing, 500, 5e5, False, pulse, station.dig, [1], 500e6, acquisition_delay_ns=0, enabled_markers=['M_SD1', 'M_SD2'])
    
    ds = do0D(fast_sweep, name='1D Coulomb Peak SD1').run()
    return ds
