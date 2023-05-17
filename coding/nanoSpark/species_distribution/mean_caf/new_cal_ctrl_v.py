import numpy as np
ctrl_v = np.empty((301,16))

for r in range(301):
    for z in range(16):
        if z == 0 or z==15:
            h = 1
        else:
            h = 2
        if r == 0:
            r1 = 0
            r3 = 1
        elif r == 300:
            r1 = 299
            r3 = 300
        else:
            r1 = r - 1
            r3 = r + 1
        ctrl_v[300-r][z] = np.pi * (r3 ** 2 - r1 ** 2) * h
np.save("./new_ctrl_v",ctrl_v)