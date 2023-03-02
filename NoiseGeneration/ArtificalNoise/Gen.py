import time

import numpy as np
import scipy.stats as ss
from matplotlib import pyplot as plt
from scipy import signal
from scipy.fftpack import next_fast_len
from scipy.signal import firwin, lfilter


def old(TargetSampleRate, df=0.01, gradient=-30, alpha=1.6, beta=0, mu=0, sigma=0.007, plot=False):
    timestart = time.time()
    np.random.seed(1)
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
    y = gradient * np.log10(freq_range) + 5  # convert to dB
    y = np.where(freq_range < 0.1, y[filterIndex], y)  # filter out low frequencies
    y = np.sqrt(10 ** (y / 10))  # convert to amplitude spectrum
    sym_y = np.zeros(len(y), dtype=float)  # make symmetric
    halfindex = int(np.ceil(len(y) / 2))  # find the middle
    sym_y[:halfindex] = y[:halfindex]  # copy the first half
    sym_y[halfindex:] = y[:halfindex - 1][::-1]  # copy the second half
    np.random.seed(1)
    GenNoise = ss.levy_stable.rvs(alpha, beta, loc=0, scale=sigma, size=len(timeSequence))  # generate noise
    fftGenNoise = np.fft.fft(GenNoise, len(GenNoise))  # fft
    fftmodified = fftGenNoise * sym_y  # apply filter
    adjustedGenNoise = np.real(np.fft.ifft(fftmodified))  # ifft
    sos = signal.butter(50, NoiseFloorCutoffFrequency, 'lowpass', analog=False, fs=fs, output="sos")  # low pass filter
    adjustedGenNoise = signal.sosfilt(sos, adjustedGenNoise)  # low pass filter
    adjustedNoiseWithNoiseFloor = adjustedGenNoise + adjustedWhiteFloorNoise
    timetaken = time.time() - timestart
    print(timetaken)
    return adjustedNoiseWithNoiseFloor


def resampled(maxf, df=0.01, gradient=-30, alpha=1.6, beta=0, mu=0, sigma=0.007, plot=False):
    dt = 1 / maxf
    maxt = round(1 / df, 2)
    timeSequence = np.arange(0, maxt + dt, dt)
    FastFFTLength = next_fast_len(len(timeSequence))
    FastFFTFrequencyRange = np.arange(0, round(FastFFTLength * df, 2), df)
    """Generate Noise floor"""
    NoiseFloorCutoffFrequency = 100
    NoiseFloorAmplitudeFactor = gradient * np.log10(NoiseFloorCutoffFrequency) - 15  # convert to dB
    NoiseFloorAmplitudeFactor = np.sqrt(10 ** (NoiseFloorAmplitudeFactor / 10))  # convert to amplitude spectrum
    WhiteFloorNoise = np.random.normal(mu, 1, FastFFTLength)  # generate noise
    fftWhiteFloorNoise = np.fft.fft(WhiteFloorNoise, len(WhiteFloorNoise))  # fft
    fftmodified = fftWhiteFloorNoise * NoiseFloorAmplitudeFactor  # apply filter
    adjustedWhiteFloorNoise = np.real(np.fft.ifft(fftmodified, len(fftmodified)))[:len(timeSequence)]  # ifft
    """ Generate Noise """
    FastFFTFrequencyRange = np.where(FastFFTFrequencyRange == 0, 1e-8, FastFFTFrequencyRange)  # avoid divide by zero
    filterIndex = np.where(FastFFTFrequencyRange > 0.1)[0][0]  # filter out low frequencies
    y = gradient * np.log10(FastFFTFrequencyRange) + 5  # convert to dB
    y = np.where(FastFFTFrequencyRange < 0.1, y[filterIndex], y)  # filter out low frequencies
    y = np.sqrt(10 ** (y / 10))  # convert to amplitude spectrum
    sym_y = np.zeros(len(y), dtype=float)  # make symmetric
    halfindex = int(np.ceil(len(y) / 2))  # find the middle
    if len(sym_y) % 2 == 0:
        sym_y[:halfindex] = y[:halfindex]  # copy the first half
        sym_y[halfindex:] = y[:halfindex][::-1]  # copy the second half
    else:
        sym_y[:halfindex] = y[:halfindex]
        sym_y[halfindex:] = y[:halfindex - 1][::-1]
    GenNoise = ss.levy_stable.rvs(alpha, beta, loc=0, scale=sigma, size=FastFFTLength)  # generate noise
    fftGenNoise = np.fft.fft(GenNoise, len(GenNoise))  # fft
    fftmodified = fftGenNoise * sym_y  # apply filter
    adjustedGenNoise = np.real(np.fft.ifft(fftmodified, len(fftmodified)))[:len(timeSequence)]  # ifft
    taps = firwin(numtaps=500, cutoff=100, width=5, pass_zero="lowpass", fs=maxf)  # low pass filter
    adjustedGenNoise = lfilter(taps, 1.0, adjustedGenNoise)  # low pass filter
    # sos = signal.butter(5000, NoiseFloorCutoffFrequency, 'lowpass', analog=False, fs=maxf, output="sos")  # low pass filter
    # adjustedGenNoise = signal.sosfilt(sos, adjustedGenNoise)  # low pass filter
    adjustedNoiseWithNoiseFloor = adjustedGenNoise + adjustedWhiteFloorNoise
    return adjustedNoiseWithNoiseFloor


if __name__ == "__main__":
    a = resampled(10000, df=8.333334e-4, gradient=-30, alpha=1.8, beta=0, mu=0, sigma=0.007,
                  plot=True)  # Compare to Real Noise
    b = old(10000, df=8.333334e-4, gradient=-30, alpha=1.8, beta=0, mu=0, sigma=0.007,
            plot=True)  # Compare to Real Noise
    dt = 1 / 10000
    maxt = 1 / 8.333334e-4
    timeSequence = np.arange(0, maxt + dt, dt)
    fig, ax = plt.subplots(3, 1)
    fig.suptitle('GeneratedNoise', fontsize=16)
    Seg_Length_mult = len(timeSequence) // 1024
    FFTLength = Seg_Length_mult * 1024
    OverlapLength = int(len(timeSequence) / 5)

    """ plot a"""
    ax[0].plot(timeSequence, a, alpha=0.5, label="Generated Noise with Noise Floor")
    ax[1].hist(a, density=True, bins="auto", histtype='stepfilled', alpha=0.5)
    ax[2].psd(a, NFFT=FFTLength, Fs=10000, alpha=0.5, noverlap=OverlapLength,
              label="Generated Noise with Noise Floor - a")
    """ plot b"""
    ax[0].plot(timeSequence, b, alpha=0.5, label="Generated Noise with Noise Floor")
    ax[1].hist(b, density=True, bins="auto", histtype='stepfilled', alpha=0.5)
    ax[2].psd(b, NFFT=FFTLength, Fs=10000, alpha=0.5, noverlap=OverlapLength,
              label="Generated Noise with Noise Floor - b")

    ax[1].set_yscale("log")
    ax[2].set_xscale("log")
    ax[2].set_xlim([1e-2, 1 / dt])
    plt.legend()
    plt.show()
