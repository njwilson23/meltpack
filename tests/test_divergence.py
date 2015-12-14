import unittest
import numpy as np
import MELTPACK.divergence

class DivergenceTests(unittest.TestCase):

    def test_mms_1(self):
        x, y = np.meshgrid(np.linspace(0, 1, 200), np.linspace(0, 1, 200))
        utest = x**2
        vtest = x*y
        htest = np.ones_like(x)

        ans = 3*x   # analytical
        div = MELTPACK.divergence.divergence(1/200, 1/200, htest, utest, vtest)
        self.assertTrue(np.mean(np.abs(ans[1:-1,1:-1]-div)) < 0.012)

    def test_mms_2(self):
        x, y = np.meshgrid(np.linspace(0, 1, 200), np.linspace(0, 1, 200))

        utest = -np.sin(5*x)
        vtest = x*y**2
        htest = np.ones_like(x)

        ans = -5*np.cos(5*x) + 2*x*y    # analytical
        div = MELTPACK.divergence.divergence(1/200, 1/200, htest, utest, vtest)
        self.assertTrue(np.mean(np.abs(ans[1:-1,1:-1]-div)) < 0.12)

if __name__ == "__main__":
    unittest.main()
