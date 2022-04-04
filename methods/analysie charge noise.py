from core_tools.sweeps.sweeps import do0D
from qcodes import MultiParameter

import numpy as np
import scipy as sp

class spectral_density(MultiParameter):
    def __init__(self, Sxx, Sxx_SD, Sxx_BG, freq):
        super().__init__('spectral_density', instrument=None, names=('Sxx','Sxx_SD', 'Sxx_BG'), shapes=((freq.size, ), (freq.size, ), (freq.size, )),
                         labels=('PDS (bg substracted)', 'Power Spectral Density (SD)', 'Power Spectral Density (background)'),
                         units=('Ev^2/Hz', 'Ev^2/Hz', 'Ev^2/Hz'),
                         setpoints=((freq,),(freq,), (freq,)),
                         setpoint_names=(('frequecy',),('frequecy',), ('frequecy',)),
                         setpoint_labels=(('Frequency',),('Frequency',), ('Frequency',)),
                         setpoint_units=(('Hz',),('Hz',), ('Hz',))
                         )
        self.data = [Sxx, Sxx_SD, Sxx_BG]
    def get_raw(self):
        return self.data

def to_PSD(x, fs, nperseg):
    n_iter = int(x.size/nperseg)
    
    freq = np.fft.rfftfreq(nperseg, d=1/fs)
    FT_coeff = np.zeros(freq.shape)
    for i in range(n_iter):
        FT_coeff += 2*np.abs(np.fft.rfft(x[i*nperseg:(i+1)*nperseg]))**2/fs/nperseg
    
    return freq, FT_coeff/n_iter

def save_PSD(obj_name, ds_1, ds_2, method = 'numpy', n_seg = 10):
        if method == 'numpy':
            freq, Sxx_SD = to_PSD( ds_1.m1(), fs = 1/ds_1.m1.x()[1], nperseg = 2**int(np.log2(int(ds_1.m1().size/n_seg))))
            freq, Sxx_BG = to_PSD( ds_2.m1(), fs = 1/ds_2.m1.x()[1], nperseg = 2**int(np.log2(int(ds_1.m1().size/n_seg))))
        if method == 'scipy':
            freq, Sxx_SD = sp.signal.welch( ds_1.m1(), fs = 1/ds_1.m1.x()[1], nperseg = 2**int(np.log2(int(ds_1.m1().size/n_seg))))
            freq, Sxx_BG = sp.signal.welch( ds_2.m1(), fs = 1/ds_2.m1.x()[1], nperseg = 2**int(np.log2(int(ds_1.m1().size/n_seg))))

        Sxx = Sxx_SD-Sxx_BG

        plt.figure()
        # plt.loglog(freq, Sxx_SD)
        # plt.loglog(freq, Sxx_BG)
        plt.loglog(freq, Sxx)
        plt.show()
        param = spectral_density(Sxx, Sxx_SD, Sxx_BG, freq)
        # return do0D(param, name=f'PSD of {obj_name}').run()
        return freq, Sxx_SD

if __name__ == '__main__':
    import matplotlib.pyplot as plt

    # # test noise funtions
    # def generate_white_noise(N, fs, P_noise):
    #   return np.random.normal(scale=np.sqrt(P_noise*fs/2), size=N)

    # fs = 1e3
    # x = generate_white_noise(100000, fs, 0.001)
    # f1, Pxx_sp = sp.signal.welch(x, fs, nperseg=10000)
    # f2, Pxx_np = to_PSD(x, fs, nperseg=10000)

    # plt.plot(f1, Pxx_sp)
    # plt.plot(f2, Pxx_np)
    # plt.show()

    # from core_tools.data.SQL.connect import set_up_local_storage
    # set_up_local_storage("xld_user", "XLDspin001", "vandersypen_data", "6dot", "XLD", "6D2S - SQ21-1-2-10-DEV-1")

    from core_tools.data.ds.data_set import load_by_uuid

    ds_1 = load_by_uuid(1648547435368974076)
    ds_2 = load_by_uuid(1648547435368974076)#this is the background

    x, y = save_PSD('SD1', ds_1, ds_2, 'numpy', n_seg=25)
    
    def one_over_f_noise(f, P_noise, alpha):
        return P_noise*f**(-alpha)
    
    def residual(param, x, data=None):
        P_noise = param['P_noise']
        alpha = param['alpha']
        
        model = one_over_f_noise(x, P_noise, alpha)

        if data is None: 
            return model
        return model-data
   
    
    from lmfit import Parameters, minimize
    
    fit_params = Parameters()
    fit_params.add('P_noise', value=1e-14)
    fit_params.add('alpha', value=1, min = 0.5, max=3)
    
    idx = np.where((x<20) & (x > 0))
    out = minimize(residual, fit_params, args=(x[idx],), kws={'data': y[idx]})
    fit = residual(out.params, x)
    plt.loglog(x, y)
    plt.loglog(x, residual(out.params, x))
    print('P_noise (eV/sqrt(Hz))',np.sqrt(out.params['P_noise']))
    print('alpha',float(out.params['alpha']))