import numpy as np
from matplotlib import pyplot as plt

poly_file = open("test_data/polynomial_set.txt").read().rstrip().split("\n")
env_file = open("test_data/lower_envelope.txt").read().rstrip().split("\n")

polynomials = []
for line in poly_file:
    data = line.split(" ")
    polynomials.append({
        "coeffs": [float(x) for x in data[0:len(data)-2]],
        "domain": [float(x) for x in data[len(data)-2:len(data)]]
    })

envelope = []
for line in env_file:
    data = line.split(" ")
    envelope.append({
        "coeffs": [float(x) for x in data[0:len(data)-2]],
        "domain": [float(x) for x in data[len(data)-2:len(data)]]
    })

def evaluate_poly(x, coeffs):
    val = 0
    for i in range(len(coeffs)):
        val += coeffs[i]*x**i
    return val

for poly in polynomials:
    x = np.linspace(poly["domain"][0], poly["domain"][1], 100)
    plt.plot(x, evaluate_poly(x, poly["coeffs"]), color="b", linestyle="dashed")

for poly in envelope:
    x = np.linspace(poly["domain"][0], poly["domain"][1], 30)
    plt.plot(x, evaluate_poly(x, poly["coeffs"]), color="r")

plt.show()