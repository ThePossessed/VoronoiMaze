import numpy as np
import matplotlib.pyplot as plt

with open("./out.txt", "rt", encoding='utf-16') as f:
    data = f.readlines()

data = [[float(j) for j in i[:-1].split(" ")] for i in data]
startEnd = data[:2]
data = data[2:]

for i in startEnd:
    plt.plot(i[0], i[1], 'go')

for i in data:
    plt.plot([i[0], i[2]], [i[1], i[3]], 'ro-')

plt.show()