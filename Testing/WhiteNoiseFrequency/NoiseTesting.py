import matplotlib.pyplot as plt
import numpy as np

Noise1 = np.random.normal(0, 0.1, 10000)
Noise2 = np.random.normal(0, 1, 10000)
FreqNoise1 = np.fft.ifft(Noise1, len(Noise1))
FreqNoise2 = np.fft.ifft(Noise2, len(Noise2))
plt.figure(1)
plt.plot(FreqNoise1)
plt.figure(2)
plt.plot(FreqNoise2)
plt.show()
