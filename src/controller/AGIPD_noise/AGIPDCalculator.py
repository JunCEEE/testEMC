#%%
import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from tqdm import tqdm
from timeit import default_timer as timer

SimExPath = '/gpfs/exfel/data/user/juncheng/simex-branch/Sources/python/'
SimExExtLib = '/gpfs/exfel/data/user/juncheng/simex-branch/lib/python3.7/'
sys.path.insert(0, SimExPath)
sys.path.insert(0, SimExExtLib)

from SimEx.Analysis.DiffractionAnalysis import DiffractionAnalysis
import SimEx.Utilities.pysingfel2dragonfly as sing2d


# %%
def gaussian(x, mu, sig):
    return 1. / (np.sqrt(2. * np.pi) * sig) * np.exp(
        -np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))


class curve_fitting:
    def __init__(self, func, xdata, ydata):
        self.__func = func
        self.__xdata = xdata
        self.__ydata = ydata
        self.__update()

    def __update(self):
        self.__popt, self.__pcov = curve_fit(self.func, self.xdata, self.ydata)
        self.__residuals = self.ydata - self.func(self.xdata, *self.popt)
        ss_res = np.sum(self.__residuals**2)
        ss_tot = np.sum((self.ydata - np.mean(self.ydata))**2)
        self.__r_squared = 1 - (ss_res / ss_tot)

    @property
    def popt(self):
        return self.__popt

    @property
    def residuals(self):
        return self.__residuals

    @property
    def r_squared(self):
        return self.__r_squared

    @property
    def func(self):
        return self.__func

    @func.setter
    def func(self, val):
        self.__func = val
        self.__update()

    @property
    def xdata(self):
        return self.__xdata

    @xdata.setter
    def xdata(self, val):
        self.__xdata = val
        self.__update()

    @property
    def ydata(self):
        return self.__ydata

    @ydata.setter
    def ydata(self, val):
        self.__ydata = val
        self.__update()

    def plotResults(self):
        xdata = self.xdata
        ydata = self.ydata
        popt = self.popt
        print('Coefficient of determination:', self.r_squared)
        for i, item in enumerate(self.popt):
            print('param {}'.format(i + 1), item)

        plt.figure()
        plt.plot(xdata, ydata, 'bo', label='data')
        plt.plot(xdata, self.func(xdata, *popt), 'g--', label='fitting')
        plt.legend()

    def predict(self, xdata):
        return self.func(xdata, *self.popt)

    def plotPredict(self, xdata):
        ydata = self.func(xdata, *self.popt)
        plt.figure()
        plt.plot(xdata, ydata)


def addBeamStop(img, stop_rad):
    """Add the beamstop in pixel radius to diffraction pattern

    :param img: Diffraction pattern
    :type img: np.2darray
    :param stop_rad: The radius of the beamstop in pixel unit
    :type stop_rad: int

    :return: Beamstop masked 2D array
    :rtype: np.2darray
    """
    stop_mask = np.ones_like(img)
    center = np.array(img.shape) // 2
    y = np.indices(img.shape)[0] - center[0]
    x = np.indices(img.shape)[1] - center[1]
    r = np.sqrt((x * x + y * y))
    stop_mask[r <= stop_rad] = 0
    masked = img * stop_mask
    return masked


def getPatternStatistic(img):
    """Get photon statistic info of a pattern

    :param img: Diffraction pattern
    :type img: np.2darray

    :return: (mean, max, min)
    :rtype: tuple
    """

    img_flat = img.ravel()
    mean_val = img_flat.mean()
    max_val = img_flat.max()
    min_val = img_flat.min()
    print('Mean: {}'.format(mean_val))
    print('Max: {}'.format(max_val))
    print('Min: {}'.format(min_val))

    return (mean_val, max_val, min_val)


def plotEstimateHist(mu, sigs_popt, n_photon=5):
    hist_xs = []
    hist_ys = []
    for i in range(n_photon):
        center = mu * i
        x_vals = np.linspace(center - 3 * linear(i, *sigs_popt),
                             center + 3 * linear(i, *sigs_popt), 120)
        y_vals = gaussian(x_vals, center, linear(i, *sigs_popt))
        hist_xs.append(x_vals)
        hist_ys.append(y_vals)
    for i in range(n_photon):
        plt.plot(hist_xs[i], hist_ys[i], label=str(i))
        plt.legend()


def getNoiseLegacy(diffr_data, mu, sigs):
    """Get the diffraction pattern with gaussian noise"""
    # zero photon noise
    n_zero = len(diffr_data[diffr_data < 1])
    n_non_zero = len(diffr_data[diffr_data >= 1])
    diffr_noise = np.zeros_like(diffr_data)
    diffr_noise[diffr_data < 1] = np.random.normal(0, sigs[0], n_zero)
    diffr_noise[diffr_data >= 1] = np.random.normal(
        mu, sigs[1], n_non_zero) * diffr_data[diffr_data >= 1]

    return diffr_noise


def getNoiseAGIPD(diffr_data, mu, sigs_popt):
    sig_arr = linear(diffr_data, *sigs_popt)
    diffr_noise = np.random.normal(diffr_data * mu, sig_arr)
    return diffr_noise


def getPopt(sigs):
    """Get the fitting parametters for predicting sigmas"""
    xdata = np.arange(len(sigs))
    ydata = sigs
    my_fitting = curve_fitting(linear, xdata, ydata)
    return my_fitting.popt


# %%
diffr_data = '/gpfs/exfel/data/user/juncheng/EMCProject/src/controller/s10/diffr.h5'
diffr_analyzer = DiffractionAnalysis(diffr_data, 1, poissonize=False)
diffr_pattern = diffr_analyzer.numpyPattern()


# %%
def linear(x, a, b):
    return a * x + b


def func_sq(x, a, b, c):
    return a * x**2 + b * x + c


# mu, sigs_popt, diffr_fn, emc_out
fwhms = np.array([49.3, 56.4, 66.4])
sigs = fwhms / 2.355
mu = 92.1  #ADU/photon
sigs_popt = getPopt(sigs)
emc_out = 'noise_data.emc'

start = timer()

diffr_analyzer = DiffractionAnalysis(diffr_data, range(1,1000), poissonize=False)
diffr_pattern = diffr_analyzer.numpyPattern()

end = timer()
print('Read pattern {} s elapsed'.format(
    end - start),flush=True)  # Time in seconds, e.g. 5.38091952400282
print('done')

# %%
diffr_noises = []
for pattern in tqdm(diffr_pattern):
    diffr_noise = getNoiseAGIPD(pattern, mu, sigs_popt)
    diffr_noises.append(diffr_noise.flatten())
# EMC sparse data
start = timer()

sPatterns = sing2d.dense_to_PatternsSOne(np.array([diffr_noise.flatten()]))
print('EMC photon shape', sPatterns.shape)
print('writing to {}'.format(emc_out))
sPatterns.write(emc_out)

end = timer()
print('{} s elapsed'.format(
    end - start))  # Time in seconds, e.g. 5.38091952400282
print('done')
# %%
