import numpy as np
import math

probabilities = []
sigma = 3.0
exp = 2.81

def genProbabilities(lower, upper):
    global probabilities
    sum = 0
    for i in range (lower, 1+upper):
        x = exp**(-math.pi*(i**2)/sigma**2)
        probabilities.append(x)
        sum += x
    probabilities[:] = [x / sum for x in probabilities]

def samplePolyGaussian(n):
    poly = []
    for i in range(0, n):
        x = np.random.choice(np.arange(-1000, 1001), p=probabilities)
        poly.append(x)
    return poly

genProbabilities(-1000, 1000)