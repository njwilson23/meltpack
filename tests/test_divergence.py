import unittest
import numpy as np
import meltpack.divergence
import matplotlib.pyplot as plt

class DivergenceTests(unittest.TestCase):

    def test_analytical_1(self):
        x, y = np.meshgrid(np.linspace(0, 1, 200), np.linspace(0, 1, 200))
        utest = x**2
        vtest = x*y
        htest = np.ones_like(x)

        ans = 3*x   # analytical
        div = meltpack.divergence.divergence(1/200, 1/200, htest, utest, vtest)
        self.assertTrue(np.mean(np.abs(ans[2:-1,2:-1]-div[1:,1:])) < 0.012)

    def test_analytical_2(self):
        x, y = np.meshgrid(np.linspace(0, 1, 200), np.linspace(0, 1, 200))

        utest = -np.sin(5*x)
        vtest = x*y**2
        htest = np.ones_like(x)

        ans = -5*np.cos(5*x) + 2*x*y    # analytical
        div = meltpack.divergence.divergence(1/200, 1/200, htest, utest, vtest)
        self.assertTrue(np.mean(np.abs(ans[2:-1,2:-1]-div[1:,1:])) < 0.12)

if __name__ == "__main__":
    unittest.main()
