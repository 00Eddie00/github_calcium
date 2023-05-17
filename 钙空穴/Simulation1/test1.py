import unittest
import numpy as np
import os

class MyTestCase(unittest.TestCase):
    def test_something(self):
        global single_area, total_area, near_triangle, num_in_triangle
        global nod, NP, gridt, npoch, NE, noe
        dir2 = "./GRID"
        i = 1
        with open(dir2 + "/gridt.dat") as file_obj:
            lines = file_obj.readlines()
            for line in lines:
                i = i + 1
        NP = i
        file_obj.close()

        i = 0
        gridt = np.zeros([2, NP], float)
        with open(dir2 + "/gridt.dat") as file_obj:
            lines = file_obj.readlines()
            for line in lines:
                line_list = line.split()
                gridt[0, i] = line_list[0]
                gridt[1, i] = line_list[1]
                i = i + 1
        file_obj.close()

        i = 0
        npoch = np.zeros(NP, int)
        with open(dir2 + "/npoch.dat") as file_obj:
            lines = file_obj.readlines()
            for line in lines:
                line_list = line.split()
                npoch[i] = line_list[0]
                i = i + 1
        file_obj.close()

        i = 1
        with open(dir2 + "/nod.dat") as file_obj:
            lines = file_obj.readlines()
            for line in lines:
                i = i + 1
        NE = i
        file_obj.close()

        i = 0
        nod = np.zeros([3, NE], int)
        with open(dir2 + "/nod.dat") as file_obj:
            lines = file_obj.readlines()
            for line in lines:
                line_list = line.split()
                nod[0, i] = line_list[0]
                nod[1, i] = line_list[1]
                nod[2, i] = line_list[2]
                i = i + 1
        file_obj.close()

        i = 0
        noe = np.zeros([3, NE], int)
        with open(dir2 + "/noe.dat") as file_obj:
            lines = file_obj.readlines()
            for line in lines:
                line_list = line.split()
                noe[0, i] = line_list[0]
                noe[1, i] = line_list[1]
                noe[2, i] = line_list[2]
                i = i + 1
        file_obj.close()

        single_area = np.zeros(NE, float)
        total_area = np.zeros(NP, float)
        cal_max = np.zeros(NP, float)
        area_matrix = np.zeros([3, 3], float)
        for i in range(0, NE):
            for j in range(0, 3):
                p = nod[j, i]
                area_matrix[j, 0] = gridt[0, p]
                area_matrix[j, 1] = gridt[1, p]
                area_matrix[j, 2] = 1.0
            area = np.linalg.det(area_matrix)
            single_area[i] = area
            for j in range(0, 3):
                p = nod[j, i]
                total_area[p] = total_area[p] + area
                cal_max[p] = cal_max[p] + 1
        nmax = max(cal_max)
        near_triangle = np.full([nmax, NP], -1)
        num_in_triangle = np.full(NP, -1)
        for i in range(0, NP):
            t = 0
            for j in range(0, NE):
                for k in range(0, 3):
                    if i == nod[k, j]:
                        near_triangle[t, i] = j
                        num_in_triangle[i] = k
                        t = t + 1

        self.assertEqual(True, False)  # add assertion here

if __name__ == '__main__':
    unittest.main()
