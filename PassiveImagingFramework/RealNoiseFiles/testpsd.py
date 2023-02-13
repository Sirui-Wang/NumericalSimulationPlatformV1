import matplotlib.pyplot as plt
import numpy as np


def plot_spectrum(s):
    f = np.fft.rfftfreq(len(s))
    return plt.loglog(f, np.abs(np.fft.rfft(s)))[0]


def noise_psd(N, psd=lambda f: 1):
    X_white = np.fft.rfft(np.random.randn(N));
    S = psd(np.fft.rfftfreq(N))
    # Normalize S
    X_shaped = X_white * S;
    plt.figure("test")
    plt.plot(X_shaped)
    return np.fft.irfft(X_shaped);


def PSDGenerator(f):
    return lambda N: noise_psd(N, f)


@PSDGenerator
def pink_noise(f):
    return 1 / np.where(f == 0, float('inf'), np.sqrt(f))


plt.style.use('dark_background')
plt.figure(figsize=(12, 8), tight_layout=True)
for G, c in zip(
        [pink_noise],
        ['hotpink']):
    plot_spectrum(G(30 * 50_000)).set(color=c, linewidth=3)
plt.legend(['pink'])
plt.suptitle("Colored Noise");
plt.ylim([1e-3, None]);
plt.show()
