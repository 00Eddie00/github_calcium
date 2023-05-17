import unittest
import numpy as np
from config.parameters import *
from math import sqrt


class MyTestCase(unittest.TestCase):
    def test_something(self):
        coefficient = np.load("config/coefficient/coefficient.npy")
        print(coefficient)

if __name__ == '__main__':
    unittest.main()
