import numpy as np
import scipy.stats as ss
import matplotlib.pyplot as plt


def normalNoise(NoiseLength, mu, sigma):
    Noise = np.random.normal(mu, sigma, NoiseLength)
    NPower = np.sum((abs(Noise)) ** 2) / len(Noise)
    Noise = Noise / np.sqrt(NPower)
    return Noise


def LevyNoise(NoiseLength, alpha, beta, mu, sigma):
    Noise = ss.levy_stable.rvs(alpha, beta, loc=mu, scale=sigma * 1 / np.sqrt(2), size=NoiseLength)
    NPower = np.sum((abs(Noise)) ** 2) / len(Noise)
    Noise = Noise / np.sqrt(NPower)
    return Noise


def RealNoise(NoiseLength):
    SensorNoise = np.random.normal(0, 0.01, NoiseLength)
    ExternalNoise = ss.levy_stable.rvs(1.75, 0, loc=0, scale=57 * 1 / np.sqrt(2), size=NoiseLength)
    AcousticNoise = ss.levy_stable.rvs(1.95, 0, loc=0, scale=10 * 1 / np.sqrt(2), size=NoiseLength)
    Noise = SensorNoise + ExternalNoise + AcousticNoise
    NPower = np.sum((abs(Noise)) ** 2) / len(Noise)
    Noise = Noise / np.sqrt(NPower)
    return Noise


def Gen(NoiseLength, NoiseType, NoiseParam=["", "", "", ""]):
    if isinstance(NoiseLength, np.ndarray):
        NoiseLength = 50 * len(NoiseLength)
    elif isinstance(NoiseLength, int):
        NoiseLength = NoiseLength
    else:
        NoiseLength = 1024
    alpha, beta, mu, sigma = NoiseParam
    if NoiseType == "Real":
        Noise = RealNoise(NoiseLength)
    elif NoiseType == "Levy":
        Noise = LevyNoise(NoiseLength, alpha, beta, mu, sigma)
    else:
        Noise = normalNoise(NoiseLength, mu, sigma)
    return Noise

# """Test"""
# dt = 0.0001
# maxt = 1000
# time = np.arange(0, maxt, dt)
# freq_range = np.arange(0,1/dt, 1/maxt)
# Noise = Gen(len(time), "Real", [1.7,0,0,50])
# plt.plot(Noise)
# plt.show()
