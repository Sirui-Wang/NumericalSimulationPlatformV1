import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as ss
from scipy import signal

np.random.seed(1107)
Sampling_Frequency = 2000  # 2kHz Sampling rate
Max_bandwidth = (0.1, 100000)  # flat response out side this range
mDuration = 600  # Measurement lasted 10 min, 600 seconds
# ss.levy_stable.parameterization="S1"

""" Distribution Parameters """
loc = 0
scale = 1
alpha = 1.9
beta = 0
mTimeSequence = np.arange(0, mDuration+1/Sampling_Frequency, 1/Sampling_Frequency)
""" Signal Generation """
normal_noise = ss.norm.rvs(loc=loc, scale=scale,size=len(mTimeSequence))
plt.figure("normal noise")
plt.plot(mTimeSequence, normal_noise)
alpha_stable_noise = ss.levy_stable.rvs(alpha, beta, loc=loc, scale=scale, size=len(mTimeSequence))
plt.figure("alpha stable noise")
plt.plot(mTimeSequence, alpha_stable_noise)

""" Compute and plot power spectral density (PSD) """
nFreqs, nPSD = signal.welch(normal_noise)
aFreqs, aPSD = signal.welch(alpha_stable_noise)

plt.figure()
# plt.semilogx(nFreqs, nPSD)
plt.semilogx(aFreqs, aPSD)
plt.show()

print("test")

NoiseLength = 4096
df = 0.1
max_f = 1000
frequency_range = np.arange(0,max_f+df, df)



def plotDist(alpha, beta, loc, scale,x_lim):
        fig, ax = plt.subplots(3,1)
        LevyStableNoise = ss.levy_stable.rvs(alpha, beta, loc=loc, scale=scale, size = NoiseLength)
        NormalNoise = ss.norm.rvs(loc=loc, scale=scale*np.sqrt(2),size = NoiseLength)
        ax_series = ax[0]
        ax_PSD = ax[1]
        ax_PDF = ax[2]
        ax_series.plot(LevyStableNoise,"r-", label="levy_stable noise time series")
        ax_series.plot(NormalNoise,"b-", label="normal noise time series")
        ax_PSD.psd(NormalNoise,NFFT=4096,Fs=2)
        # Levy_PSD = np.fft.fftshift(np.fft.fft(LevyStableNoise))
        # Levy_PSD = Levy_PSD[NoiseLength//2:] # only look at positive frequencies.  remember // is just an integer divide
        # Normal_PSD = np.fft.fftshift(np.fft.fft(NormalNoise))
        # Normal_PSD = Normal_PSD[NoiseLength//2:] # only look at positive frequencies.  remember // is just an integer divide
        # ax_PSD.plot(Levy_PSD,"r-", label="levy_stable noise time series")
        # ax_PSD.plot(Normal_PSD,"b-", label="normal noise time series")



def main():
        np.random.seed(1)
        alpha = 1
        beta = 0
        loc = 2
        scale = 1
        x_lim =(-10+loc,10+loc)
        plotDist(alpha, beta, loc, scale, x_lim)

main()
plt.show()

# alpha = 1  # characteristic exponents that discribes the tail of the distribution. （0，2] -- 0<alpha<=2
# beta = -1 # skewness that defines if the distribution is left skewed (beta <0) or right skewed (beta >0) [-1,1]
# loc = 0 # similar to mean in normal distribution
# scale = 1 # similar to std in normal distribution
# # x = np.linspace(ss.levy_stable.ppf(-0.99, alpha, beta,loc=loc, scale=scale),
# #                 ss.levy_stable.ppf(0.99, alpha, beta,loc=loc, scale=scale), 100)
# # x_norm = np.linspace(ss.norm.ppf(0.001, loc=loc, scale=scale*np.sqrt(2)), ss.norm.ppf(0.999,loc=loc, scale=scale*np.sqrt(2)),100)
# x = np.arange(-4,4.01,0.1)
# ASSeries = ss.levy_stable.rvs(alpha, beta, loc=loc, scale=scale, size = 1000)
# NSeries = ss.norm.rvs(loc=loc, scale=scale*np.sqrt(2), size=1000)
# fig, ax = plt.subplots(2,1)
# ax[0].plot(ASSeries,"r-", label="levy_stable noise time series")
# ax[0].plot(NSeries,"b-", label="normal noise time series")
# ax[1].hist(ASSeries, density=True, bins="auto", histtype='stepfilled', alpha=0.5)
# ax[1].set_xlim([x[0], x[-1]])
# ax[1].plot(x, ss.norm.pdf(x, loc=loc, scale=scale*np.sqrt(2)),
#         'b-', lw=5, alpha=0.6, label='normal pdf')
# ax[1].plot(x, ss.levy_stable.pdf(x, alpha, beta, loc=loc, scale=scale),
#         'r-', lw=1, alpha=1, label='levy_stable pdf')
# ax[0].legend(loc='best', frameon=False)
# ax[1].legend(loc='best', frameon=False)
plt.show()
