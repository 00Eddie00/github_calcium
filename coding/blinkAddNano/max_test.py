import unittest
import numpy as np


class MyTestCase(unittest.TestCase):
    def test_something(self):
        ampl = np.asarray([1, 2, 3, 4, 5, 6, 7, 8, 9])
        peak_value = np.max(ampl)  # 峰值，即某空间点在时间尺度上的最大荧光强度
        print("peak_value:", peak_value)
        t_rise_index = np.argmax(ampl)  # 峰化时间对应ampl中的索引
        print("t_rise_index:", t_rise_index)
        temp = np.maximum(ampl - peak_value / 2, peak_value / 2 - ampl)
        print("temp:", temp)
        t50_index = np.argmin(temp[t_rise_index:])  # 从最大衰减到50%时索引
        print("t50_index:", t50_index)
        fdhm_index = np.argmin(temp[:t_rise_index])  # 从50%增加到最大时索引
        print("fdhm_index:", fdhm_index)


if __name__ == '__main__':
    unittest.main()
