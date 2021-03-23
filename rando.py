import random
import numpy as np
vals = []

for i in range(100000):
    vals.append(np.round(np.random.normal(0, scale=1/3)))

val = set(vals)

for x in val:
    print(x," ", vals.count(x))

print("Total ", len(val))