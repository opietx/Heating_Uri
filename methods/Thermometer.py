from good_morning.calibrations.ultility import get_target, readout_convertor

from core_tools.utility.variable_mgr.var_mgr import variable_mgr
from core_tools.sweeps.pulse_lib_wrappers.PSB_exp import run_qubit_exp
from core_tools.sweeps.sweeps import scan_generic
from core_tools.job_mgnt.job_mgmt import job_wrapper

from dev_V2.six_qubit_QC_v2.system import six_dot_sample
from dev_V2.Elzerman_2_qubits_clean.TRIG import mk_TRIG
from dev_V2.six_qubit_QC_v2.VAR import variables

import pulse_lib.segments.utility.looping as lp
import matplotlib.pyplot as plt

import qcodes as qc
import scipy as sp
import numpy as np
from core_tools.utility.variable_mgr.var_mgr import variable_mgr

from core_tools.GUI.keysight_videomaps.data_getter.scan_generator_Keysight import  construct_1D_scan_fast, construct_2D_scan_fast
from core_tools.sweeps.sweeps import do0D, do1D


from qcodes import Parameter
import scipy



# @job_wrapper
def Addition_line1(sweep_range=25, plot=False):
    var_mgr = variable_mgr()


    fast_sweep_1D = construct_1D_scan_fast('vP1', sweep_range, 200,5e5, True, pulse, station.dig, [1], 500e6, pulse_gates={'vP1': var_mgr.PSB12_P1+3}, acquisition_delay_ns=50e3, enabled_markers=['M_SD1', 'M_SD2'])

    init= station.MW_source.frequency()
    ds = do1D(station.MW_source.frequency,init ,init+10,100,1e-3, fast_sweep_1D, name='1D addition line',reset_param=True).run()
    return ds
def Addition_line2(sweep_range=25, plot=False):
    var_mgr = variable_mgr()


    fast_sweep_1D = construct_1D_scan_fast('vP2', sweep_range, 200,5e5, True, pulse, station.dig, [1], 500e6, pulse_gates={'vP2': var_mgr.PSB12_P2}, acquisition_delay_ns=50e3, enabled_markers=['M_SD1', 'M_SD2'])

    init= station.MW_source.frequency()
    ds = do1D(station.MW_source.frequency,init ,init+10,100,1e-3, fast_sweep_1D, name='1D addition line',reset_param=True).run()
    return ds



def FermiLinear(x, a, b, cc, A, T, l=1.16):
    """ Fermi distribution with linear function added
    Arguments:
        x (numpy array): independent variable
        a, b (float): coefficients of linear part
        cc (float): center of Fermi distribution
        A (float): amplitude of Fermi distribution
        T (float): temperature Fermi distribution in Kelvin
        l (float): leverarm divided by kb
    The default value of the leverarm is
        (100 ueV/mV)/kb = (100*1e-6*scipy.constants.eV )/kb = 1.16.
    For this value the input variable x should be in mV, the
    temperature T in K. We input the leverarm divided by kb for numerical stability.
    Returns:
        y (numpy array): value of the function
    .. math::
        y = a*x + b + A*(1/ (1+\exp( l* (x-cc)/(T) ) ) )
    """
    y = a * x + b + A * 1. / (1 + np.exp(l * (x - cc) / (T)))
    return y


def linear_function(x,a,b):
    return x*a + b


def _estimate_fermi_model_center_amplitude(x_data, y_data_linearized):
    
    """ Estimates the following properties of a charge addition line; the center location
        of the addition line. The amplitude step size caused by the addition line.
        
        Borrowed from QTT
        
    Args:
        x_data (1D array): The independent data.
        y_data_linearized (1D array): The dependent data with linear estimate subtracted.
    Returns:
        xdata_center_est (float): Estimate of x-data value at the center.
        amplitude_step (float): Estimate of the amplitude of the step.
    """
    sigma = x_data.size / 250
    y_derivative_filtered = scipy.ndimage.gaussian_filter(y_data_linearized, sigma, order=1)

    # assume step is steeper than overall slope
    estimated_index = np.argmax(np.abs(y_derivative_filtered))
    center_index = int(x_data.size / 2)

    # prevent guess to be at the edges
    if estimated_index < 0.01 * x_data.size or estimated_index > 0.99 * x_data.size:
        estimated_center_xdata = np.mean(x_data)
    else:
        estimated_center_xdata = x_data[estimated_index]

    split_offset = int(np.floor(x_data.size / 10))
    mean_right = np.mean(y_data_linearized[(center_index + split_offset):])
    mean_left = np.mean(y_data_linearized[:(center_index - split_offset)])
    amplitude_step = -(mean_right - mean_left)

    # if np.sign(-y_derivative_filtered[estimated_index]) != np.sign(amplitude_step):
        # warnings.warn('step size might be incorrect')

    return estimated_center_xdata, amplitude_step


def initFermiLinear(x_data, y_data):
    """ Initialization of fitting a FermiLinear function.
    First the linear part is estimated, then the Fermi part of the function.
    
    This funciton in 99% obtained form QTT
    
    Args:
        x_data (array): data for independent variable
        y_data (array): dependent variable
    Returns:
        linear_part (array)
        fermi_part (array)
    """
    xdata = np.array(x_data)
    ydata = np.array(y_data)
    n = xdata.size
    nx = int(np.ceil(n / 5))

    if nx < 4:
        p1, _ = scipy.optimize.curve_fit(linear_function, np.array(xdata[0:100]),
                                         np.array(ydata[0:100]))

        a = p1[0]
        b = p1[1]
        linear_part = [a, b]
        ylin = linear_function(xdata, linear_part[0], linear_part[1])
        cc = np.mean(xdata)
        A = 0
        T = np.std(xdata) / 10
        fermi_part = [cc, A, T]
    else:
        # guess initial linear part
        mx = np.mean(xdata)
        my = np.mean(ydata)
        dx = np.hstack((np.diff(xdata[0:nx]), np.diff(xdata[-nx:])))
        dx = np.mean(dx)
        dd = np.hstack((np.diff(ydata[0:nx]), np.diff(ydata[-nx:])))
        dd = np.convolve(dd, np.array([1., 1, 1]) / 3)  # smooth
        if dd.size > 15:
            dd = np.array(sorted(dd))
            w = int(dd.size / 10)
            a = np.mean(dd[w:-w]) / dx
        else:
            a = np.mean(dd) / dx
        b = my - a * mx
        linear_part = [a, b]
        ylin = linear_function(xdata, *linear_part)

        # subtract linear part
        yr = ydata - ylin

        cc, A = _estimate_fermi_model_center_amplitude(xdata, yr)

        T = np.std(xdata) / 100
        linear_part[1] = linear_part[1] - A / 2  # correction
        fermi_part = [cc, A, T]

        yr = ydata - linear_function(xdata, *linear_part)
    return linear_part, fermi_part


def Fermi_center(x_data, y_data, plot=False):
    '''
    Finds the center of a Fermi distribution via fitting. First makes a guess and then fits to return only the center
    Arguments:
        x_data (array) = x axis of measurement, the AWG-swing voltage setpoints
        y_data (array) = y-axis of measurement, the measured voltage on digitizer
    Returns:
        center value of the fermi distribution
    '''
    # plt.figure()
    # plt.plot(x_data,y_data, '.')
    
    '''Initial guess'''
    pp_linera,pp_femri=initFermiLinear(x_data,y_data)
    initial_parameters=[pp_linera[0],pp_linera[1],pp_femri[0],pp_femri[1],pp_femri[2]]
    # y_guess = FermiLinear(x_data,*initial_parameters)
    # plt.plot(x_data,y_guess, '--')

    '''fit'''
    fitting_results,_ = scipy.optimize.curve_fit(FermiLinear, x_data, y_data, p0=initial_parameters)
    
    y_fit = FermiLinear(x_data,*fitting_results)
    if plot == True:
        plt.figure()
        plt.plot(x_data,y_data,'o')
        plt.plot(x_data,y_fit,linewidth=5, alpha=0.5)
        plt.show()
    return fitting_results

#%%


if __name__ == "__main__":
    # ds = Addition_line1(25)
    # ds = Addition_line2(25)
    
    ds_avg = ds['ch1']
    a= Fermi_center(ds_avg.x()[1:],ds_avg.y()[1:], plot=True)
    
    

