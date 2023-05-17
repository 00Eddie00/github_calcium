import matplotlib.pyplot as plt
from optical_blurring.process_concentrationv2 import *

def fluo_contourf(result_set, species, c_val):
    # 准备数据
    r = np.linspace(-300, 300, 601)  # r轴坐标
    z = np.linspace(-7.5, 7.5, 16)  # z轴坐标
    # 处理浓度数据
    original_data = np.loadtxt(f"{result_set}/{species}/{species}00010000.csv")
    con_matrix = process_concentration(original_data)[300,:,:].T
    # 绘制
    plt.contourf(r, z, con_matrix, levels=1000)
    plt.xlabel('r')  # 设置x轴标签
    plt.ylabel('z')  # 设置y轴标签
    plt.colorbar()
    # 显示图表
    plt.show()

fluo_contourf("../../Result/TS100000_更改迭代BUG后","CaF",0)