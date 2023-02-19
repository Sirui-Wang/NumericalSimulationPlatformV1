import numpy as np
from matplotlib import animation
from matplotlib import pyplot as plt

fig = plt.figure()
ax = plt.axes(xlim=(0, 2), ylim=(-2, 2))

line, = ax.plot([], [])


def animate(i):
    x = np.linspace(0, 2, 1000)
    y = np.sin(2 * np.pi * (x - 0.01 * i))
    y -= y % 0.3
    line.set_data(x, y)
    return line,


anim = animation.FuncAnimation(fig, animate,
                               frames=100, interval=100, blit=True)

plt.show()
