import matplotlib.pyplot as plt
from math import *
import numpy as np

measurements = np.array([
    1.00,
    1.07,
    1.10,
    1.10,
    1.16,
    1.20,
    1.22,
    1.29,
    1.29,
    1.31,
    1.34,
    1.39,
    1.41,
    1.48,
    1.50,
    1.53,
    1.74,
    1.72,
    1.72,
    1.81
])

formula = np.array([
    1.00,
    1.10,
    1.14,
    1.14,
    1.25,
    1.29,
    1.31,
    1.38,
    1.36,
    1.39,
    1.58,
    1.62,
    1.71,
    1.75,
    1.77,
    1.91,
    2.11,
    1.97,
    2.00,
    2.22
])

array_to_use = formula / measurements

n = 14

hot = np.array([255, 32, 0])
cold = np.array([0, 32, 255])


def getcolor(v: float):
    l = (v - min(array_to_use)) / (max(array_to_use) - min(array_to_use))
    d = l * hot + (1.0 - l) * cold
    c = int(int(d[0]) * 16**4 + int(d[1]) * 16**2 + int(d[2]))
    s = hex(c)[2:]
    while len(s) < 6:
        s = "0" + s
    return s

# FOR PRINTING A TABLE


# k = 0

# print(" & ", end="")
# for i in range(n):
#     print("{}".format(i + 1), end="")
#     if i != n - 1:
#         print(" & ", end="")

# print("\\\\")
# print("\\hline")

# for i in range(0, n):
#     print("{} & ".format(i + 1), end="")
#     for j in range(0, i):
#         print(" & ", end="")
#     for j in range(i, n):
#         print("\\textcolor[HTML]{", getcolor(array_to_use[k]),  "}", end="")
#         print("{{{:.2f}}}".format(array_to_use[k], 2), end="")
#         if j != n - 1:
#             print(" & ", end="")
#         k += 1
#     print(" \\\\")

# FOR PRINTING AN ARRAY


for i in range(0, len(array_to_use)):
    print("\\textcolor[HTML]{", getcolor(array_to_use[i]),  "}", end="")
    print("{{{:.2f}}}".format(array_to_use[i], 2))
