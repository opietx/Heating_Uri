from core_tools.sweeps.sweeps import do0D, do1D, do2D, scan_generic
import matplotlib.pyplot as plt
import matplotlib
import time

def get_slope_and_get_meas_points(ds, background_subtraction= True):
    """ Allows the user to click three points. This method will calculate the (positive)
        slope between point 1 and 2 and defines the coordinates of point 3.
        Args:
            ds (dataset) : dataset with charge stability diagram and gate voltage in mV
            background_subtraction (bool) : include background subtraction as clickable point
        Returns:
            slope (double) : the slope between point 1 and point 2 (always positive)
            noise_meas (double) : the x coordinate of the noise measurement
            bg_meas (double) : the x coordinate of the background measurement
    """
    fig = plt.figure()
    plt.plot(ds.m1.x(),ds.m1.y())
    plt.xlabel(f'{ds.m1.x.label} ({ds.m1.x.unit})')
    plt.ylabel(f'{ds.m1.y.label} ({ds.m1.y.unit})')

    print('''Instructions ::\n\n
        1) click two points along the line that you want to fit.\n
        2) Then click the point at which you want to measure on the slope.\n 
        3) Take the fourth point on Coulomb blockade as calibration.''')
    
    coords = []
    
    def onclick(event):
        ix, iy = event.xdata, event.ydata
        coords.append((ix, iy))
        
        plt.scatter(ix, iy, s=80, c='r')
        fig.canvas.draw()
        
        if len(coords) == 2:
            print(f'slope measured : {(coords[1][1]-coords[0][1])/(coords[1][0]-coords[0][0])} [A/V]')
        if len(coords) == 3:
            print(f'measure noise at {ix} [mV]')
        if len(coords) == 4:
            print(f'measure background at {ix} [mV]')
            fig.canvas.mpl_disconnect(cid)
            plt.close(fig)
            
    cid = fig.canvas.mpl_connect('button_press_event', onclick)
    plt.show()
    
    while len(coords) != 4:
        plt.pause(1)

    slope = (coords[1][1]-coords[0][1])/(coords[1][0]-coords[0][0])*1000 #unit in A/V [gate in mV]
    to_meas = coords[2][0]
    bg_meas = coords[3][0]

    return slope, to_meas, bg_meas

def measure_slope(gate, v_start, v_stop, m_instr, lever_arm, background_subtraction = False):
    '''
    Measure the noise on a keithley DMM6500
    
    Args:
        gate (paramter) : parameter of qcodes
        v_start (double) : start voltage (mV)
        v_stop (double) : stop voltage (mV)
        m_instr (QCoDeS Instrument) : assumsed keithley6500 driver in core tools
        lever_arm (double) : lever arm of the gate to the dot (eV/V)

    Returns:
        time_trace_noise, time_trace_background (dataset) : data sets of the noise under the different conditions.
    '''
    # measure coulomb peak
    original_val = gate()
    ds = do1D(gate, v_start, v_stop, 200, 0.01, m_instr.amplitude, name = f'1D gate scan of {gate.name}').run()

    # ask user to fit the slope
    slope, to_meas, bg_meas = get_slope_and_get_meas_points(ds, background_subtraction=background_subtraction)

    # set gate to the noisy point and measure a time trace
    gate(to_meas)
    time.sleep(1)
    m_instr.prepare_buffer(600, 1200000)
    trace_obj = m_instr.TimeTrace
    trace_obj.multiplier = lever_arm/slope
    trace_obj.labels = ('Energy',)
    trace_obj.units = ('eV',)
    time_trace_noise = do0D(trace_obj, name = f'time trace {gate.name}').run()
    
    # set gate to the background point and measure a time trace
    gate(bg_meas)
    time.sleep(1)
    m_instr.prepare_buffer(600, 1200000)
    trace_obj = m_instr.TimeTrace
    trace_obj.multiplier = lever_arm/slope
    trace_obj.labels = ('Energy',)
    trace_obj.units = ('eV',)
    time_trace_background = do0D(trace_obj, name = f'time trace {gate.name} bg_signal').put()

    # wait 1 second to allow the gate voltage to stabilize
    gate(original_val)
    return time_trace_noise, time_trace_background

#%%
if __name__ == "__main__":    
    gate = station.gates.vSD1_P
    
    v_start = gate()-7.5
    v_stop =  gate()+7.5
    inst= [station.Idc_1]
    ds = do1D(gate, v_start, v_stop, 200, 0.01, *inst,reset_param=True,name = f'1D gate scan of {gate.name}').run()
    # measure_slope(gate, v_start, v_stop, station.keithley_3, 0.013)

