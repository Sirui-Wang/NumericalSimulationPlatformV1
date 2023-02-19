import matplotlib.pyplot as plt
import numpy
import numpy as np
import pandas as pd
import scipy.stats as ss
from scipy import signal


def resampled(TargetSampleRate, df=0.01, gradient=-30, alpha=1.6, beta=0, mu=0, sigma=0.007, plot=False):
    fs = TargetSampleRate
    dt = 1 / TargetSampleRate
    maxt = 1 / df
    timeSequence = np.arange(0, maxt + dt, dt)
    freq_range = np.arange(0, TargetSampleRate + df, df)
    """Generate Noise floor"""
    NoiseFloorCutoffFrequency = 100
    NoiseFloorAmplitudeFactor = gradient * np.log10(NoiseFloorCutoffFrequency) - 15  # convert to dB
    NoiseFloorAmplitudeFactor = np.sqrt(10 ** (NoiseFloorAmplitudeFactor / 10))  # convert to amplitude spectrum
    WhiteFloorNoise = np.random.normal(mu, 1, len(timeSequence))  # generate noise
    fftWhiteFloorNoise = np.fft.fft(WhiteFloorNoise, len(WhiteFloorNoise))  # fft
    fftmodified = fftWhiteFloorNoise * NoiseFloorAmplitudeFactor  # apply filter
    adjustedWhiteFloorNoise = np.real(np.fft.ifft(fftmodified))  # ifft
    """ Generate Noise """
    freq_range = np.where(freq_range == 0, 1e-8, freq_range)  # avoid divide by zero
    filterIndex = np.where(freq_range > 0.1)[0][0]  # filter out low frequencies
    y = gradient * np.log10(freq_range)  # convert to dB
    y = np.where(freq_range < 0.1, y[filterIndex], y)  # filter out low frequencies
    y = np.sqrt(10 ** (y / 10))  # convert to amplitude spectrum
    sym_y = np.zeros(len(y), dtype=float)  # make symmetric
    halfindex = int(np.ceil(len(y) / 2))  # find the middle
    sym_y[:halfindex] = y[:halfindex]  # copy the first half
    sym_y[halfindex:] = y[:halfindex - 1][::-1]  # copy the second half
    GenNoise = ss.levy_stable.rvs(alpha, beta, loc=0, scale=sigma, size=len(timeSequence))  # generate noise
    fftGenNoise = np.fft.fft(GenNoise, len(GenNoise))  # fft
    fftmodified = fftGenNoise * sym_y  # apply filter
    adjustedGenNoise = np.real(np.fft.ifft(fftmodified))  # ifft
    # N, Wn = signal.buttord(0.2,0.3,3,200,analog=False,fs=fs)
    # print(N, Wn)
    sos = signal.butter(50, NoiseFloorCutoffFrequency, 'lowpass', analog=False, fs=fs, output="sos")  # low pass filter
    adjustedGenNoise = signal.sosfilt(sos, adjustedGenNoise)  # low pass filter
    adjustedNoiseWithNoiseFloor = adjustedGenNoise + adjustedWhiteFloorNoise
    if plot:
        fig, ax = plt.subplots(3, 1)
        fig.suptitle('GeneratedNoise', fontsize=16)
        """ plot generated Noise """
        Seg_Length_mult = len(freq_range) // 1024
        FFTLength = Seg_Length_mult * 1024
        OverlapLength = int(len(freq_range) / 5)
        # ax[0].plot(timeSequence, adjustedGenNoise / max(abs(adjustedGenNoise)), label="Generated Noise")
        # ax[1].hist(adjustedGenNoise, density=True, bins="auto", histtype='stepfilled')
        # ax[2].psd(adjustedGenNoise, NFFT=FFTLength, Fs=fs, window=numpy.blackman(FFTLength), noverlap=OverlapLength, label="Generated Noise")
        """Compare to Generated Noise with Noise Floor"""
        ax[0].plot(timeSequence, adjustedNoiseWithNoiseFloor, alpha=0.5, label="Generated Noise with Noise Floor")
        ax[1].hist(adjustedNoiseWithNoiseFloor, density=True, bins="auto", histtype='stepfilled', alpha=0.5)
        ax[1].set_yscale("log")
        Seg_Length_mult = len(freq_range) // 1024
        FFTLength = Seg_Length_mult * 1024
        OverlapLength = int(len(freq_range) / 5)
        ax[2].psd(adjustedNoiseWithNoiseFloor, NFFT=FFTLength, Fs=fs, alpha=0.5, window=numpy.blackman(FFTLength),
                  noverlap=OverlapLength, label="Generated Noise with Noise Floor")
        ax[2].set_xscale("log")
        ax[2].set_xlim([1e-2, 1 / dt])
        """Compare to Real Noise"""
        NoiseArray = pd.read_csv("../../Data/RealNoiseSequence.csv", delimiter=",", engine="c").to_numpy().reshape(-1)[
                     :12000001]
        ax[0].plot(timeSequence, NoiseArray, alpha=0.5, label="Real Noise")
        ax[1].hist(NoiseArray, density=True, bins="auto", histtype='stepfilled', alpha=0.5)
        ax[2].psd(NoiseArray, NFFT=FFTLength, Fs=fs, window=numpy.blackman(FFTLength), alpha=0.5,
                  noverlap=OverlapLength, label="Real Noise")
        plt.legend()
        plt.show()
    return adjustedNoiseWithNoiseFloor


if __name__ == "__main__":
    resampled(10000, df=8.333334e-4, gradient=-30, alpha=1.6, beta=0, mu=0, sigma=0.007,
              plot=True)  # Compare to Real Noise

#
# import matplotlib.pyplot as plt
# import numpy
# import numpy as np
# import pandas as pd
# import scipy.stats as ss
# from scipy.stats import levy_stable
# from scipy import signal
# from tkinter import *
# from tkinter import filedialog
# import os
#
#
# def resampled(TargetSampleRate, df=0.01, gradient=-30, alpha=1.6, beta=0, mu=0, sigma=0.007, seed=False):
#     if seed:
#         np.random.seed(1107)
#     fs = TargetSampleRate
#     dt = 1/TargetSampleRate
#     maxt = 1/df
#     timeSequence = np.arange(0, maxt+dt, dt)
#     freq_range = np.arange(0, TargetSampleRate+df, df)
#     """Generate Noise floor"""
#     NoiseFloorCutoffFrequency = 100
#     NoiseFloorAmplitudeFactor = gradient * np.log10(NoiseFloorCutoffFrequency)-20 # convert to dB
#     NoiseFloorAmplitudeFactor = np.sqrt(10 ** (NoiseFloorAmplitudeFactor / 10))  # convert to amplitude spectrum
#     WhiteFloorNoise = np.random.normal(mu, 1, len(timeSequence))  # generate noise
#     fftWhiteFloorNoise = np.fft.fft(WhiteFloorNoise, len(WhiteFloorNoise))  # fft
#     fftWhiteNoiseModified = fftWhiteFloorNoise * NoiseFloorAmplitudeFactor  # apply filter
#     adjustedWhiteFloorNoise = np.real(np.fft.ifft(fftWhiteNoiseModified))  # ifft
#     """ Generate Noise """
#     freq_range = np.where(freq_range == 0, 1e-8, freq_range) # avoid divide by zero
#     filterIndex = np.where(freq_range > 0.1)[0][0] # filter out low frequencies
#     y = gradient * np.log10(freq_range) - 40  # convert to dB
#     y = np.where(freq_range < 0.1, y[filterIndex], y) # filter out low frequencies
#     y = np.sqrt(10 ** (y / 10))  # convert to amplitude spectrum
#     sym_y = np.zeros(len(y), dtype=float)  # make symmetric
#     halfindex = int(np.ceil(len(y) / 2))  # find the middle
#     sym_y[:halfindex] = y[:halfindex]  # copy the first half
#     sym_y[halfindex:] = y[:halfindex - 1][::-1]  # copy the second half
#     GenNoise = ss.levy_stable.rvs(alpha, beta, loc=mu, scale=sigma, size=len(timeSequence))  # generate noise
#     fftGenNoise = np.fft.fft(GenNoise, len(GenNoise))  # fft
#     fftmodified = fftGenNoise * sym_y # apply filter
#     adjustedGenNoise = np.real(np.fft.ifft(fftmodified))  # ifft
#     sos = signal.butter(10, NoiseFloorCutoffFrequency, 'lowpass', analog=False, fs=fs, output="sos")  # low pass filter
#     adjustedGenNoise = signal.sosfilt(sos, adjustedGenNoise)  # low pass filter
#     adjustedGenNoise += adjustedWhiteFloorNoise
#     fig, ax = plt.subplots(3, 1)
#     fig.suptitle('GeneratedNoise', fontsize=16)
#     """ plot generated Noise """
#     Seg_Length_mult = len(freq_range) // 1024
#     """Compare to Generated Noise with Noise Floor"""
#     ax[0].plot(timeSequence, adjustedGenNoise, alpha=0.5, label="Generated Noise with Noise Floor")
#     ax[1].hist(adjustedGenNoise, density=True, bins="auto", histtype='stepfilled', alpha=0.5)
#     Seg_Length_mult = len(freq_range) // 1024
#     FFTLength = Seg_Length_mult * 1024
#     OverlapLength = int(len(freq_range) / 5)
#     ax[2].psd(adjustedGenNoise, NFFT=FFTLength, Fs=fs, alpha = 0.5, window=numpy.blackman(FFTLength), noverlap=OverlapLength, label="Generated Noise with Noise Floor")
#     ax[2].set_xscale("log")
#     ax[2].set_xlim([1e-2, 1 / dt])
#     """Compare to Real Noise"""
#     NoiseArray = pd.read_csv("../../Data/RealNoiseSequence.csv", delimiter=",", engine="c").to_numpy().reshape(-1)[:12000001]
#     ax[0].plot(timeSequence, NoiseArray, alpha=0.5, label="Real Noise")
#     ax[1].hist(NoiseArray, density=True, bins="auto", histtype='stepfilled', alpha=0.5)
#     ax[2].psd(NoiseArray, NFFT=FFTLength, Fs=fs, window=numpy.blackman(FFTLength), alpha = 0.5,noverlap=OverlapLength,label="Real Noise")
#     plt.legend()
#     plt.show()
#     return adjustedGenNoise
#
#
#
# resampled(10000, df=8.333334e-4, gradient=-30, alpha=1.6, beta=0, mu=0, sigma=0.007, seed=True) # Compare to Real Noise
