import numpy as np

parent_path = "E:\\Code\\Python\\PycharmProjects\\blinkAddNano\\blinkgrid\\"

point_och = np.loadtxt(parent_path + "npoch.dat", dtype=np.int64)

for i in range(len(point_och)):
    # 1 壁面  2 入流  4 出流
    if(point_och[i]  == 4):
        point_och[i] = 2
    # 将入流边界与壁面边界改为出流边界
    elif(point_och[i] == 2 or point_och[i] == 1):
        point_och[i] = 4

np.savetxt(f"{parent_path}npoch_nano.dat", point_och, "%d")