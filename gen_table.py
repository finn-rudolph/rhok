import matplotlib.pyplot as plt
from math import *
import numpy as np

y = []
n = 14

hot = np.array([255, 32, 0])
cold = np.array([0, 32, 255])


def getcolor(v: float):
    l = (v - min(y)) / (max(y) - min(y))
    d = l * hot + (1.0 - l) * cold
    c = int(int(d[0]) * 16**4 + int(d[1]) * 16**2 + int(d[2]))
    s = hex(c)[2:]
    while len(s) < 6:
        s = "0" + s
    return s


k = 0

print(" & ", end="")
for i in range(n):
    print("{}".format(i + 1), end="")
    if i != n - 1:
        print(" & ", end="")

print("\\\\")
print("\\hline")

for i in range(0, n):
    print("{} & ".format(i + 1), end="")
    for j in range(0, i):
        print(" & ", end="")
    for j in range(i, n):
        print("\\textcolor[HTML]{", getcolor(y[k]),  "}", end="")
        print("{{{:.2f}}}".format(y[k], 2), end="")
        if j != n - 1:
            print(" & ", end="")
        k += 1
    print(" \\\\")
