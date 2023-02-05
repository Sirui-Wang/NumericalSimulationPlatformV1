from scipy.stats import levy_stable
import matplotlib.pyplot as plt

points = 1000
jennys_constant = 8675309
alpha, beta = 1.98, -0.5

draw = levy_stable.rvs(alpha, beta, size=points, random_state=jennys_constant)
plt.plot(draw)
plt.show()
