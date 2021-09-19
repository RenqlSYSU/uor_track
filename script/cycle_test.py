import matplotlib.pyplot as plt
from matplotlib import animation
import numpy as np

a = np.random.rand(10,10)

fig, ax=plt.subplots()
container = []

for i in range(a.shape[1]):
    line,  = ax.plot(a[:,i])
    title = ax.text(0.5,1.05,"Title {}".format(i), 
                    size=plt.rcParams["axes.titlesize"],
                    ha="center", transform=ax.transAxes, )
    container.append([line, title])

ani = animation.ArtistAnimation(fig, container, interval=200, blit=False)

plt.show()
ani.save('./2.mp4', fps=30)

