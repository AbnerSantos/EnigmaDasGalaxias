import matplotlib.pyplot as plt
import numpy as np

n = int(input())

galaxies_x = np.zeros(n)
galaxies_y = np.zeros(n)

for i in range(n):
    point = input().split()
    galaxies_x[i] = int(point[0])
    galaxies_y[i] = int(point[1])

plt.scatter(galaxies_x, galaxies_y, c='b')
plt.show()
