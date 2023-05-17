import numpy as np
from config.parameters import R0, R1, NR


def generate_spark_mesh():
    """
    生成网格信息，并写入config目录下spark_mesh.csv文件中
    :return:
    :rtype:
    """
    spark_mesh = np.linspace(R0, R1, NR)
    np.savetxt("../config/spark_mesh.csv", spark_mesh, delimiter=",")


def load_spark_mesh():
    """
    加载网格信息，并返回
    :return: 网格信息
    :rtype: ndarray NR
    """
    spark_mesh = np.loadtxt("../config/spark_mesh.csv", delimiter=",")
    return spark_mesh

if __name__ == '__main__':
    generate_spark_mesh()