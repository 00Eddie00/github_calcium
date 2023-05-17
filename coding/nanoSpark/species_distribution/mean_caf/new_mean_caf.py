import numpy as np
import os

path = "../../Result/TS70000_kRyR=311999711.0000832_09-23-11-38-22-new_grid/convolved_result"
dir_list = os.listdir(path)
mean_caf = np.empty(len(dir_list))

ctrl_v = np.load("./new_ctrl_v.npy")
total_v = sum(sum(ctrl_v))

count = 0
for i in dir_list:
    print(count)
    temp_caf = np.loadtxt(f"{path}/{i}")[0:301]
    mean_caf[count] = sum(sum(temp_caf * ctrl_v)) / total_v
    count += 1

np.savetxt("./mean_caf.csv", mean_caf)