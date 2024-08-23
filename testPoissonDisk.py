import numpy as np
import matplotlib.pyplot as plt

with open("./out.txt", "rt", encoding='utf-16') as f:
    data = f.readlines()

data = [[float(j) for j in i[:-1].split(" ")] for i in data]

for i in data:
    plt.plot(i[0], i[1], 'go')

plt.show()