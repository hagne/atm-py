import unittest

import numpy as np
import pandas as pd

from atmPy.radiation.retrievals.langley_calibration import fit_langley


class FitLangleyTest(unittest.TestCase):
    def test_robust_weighted_fit(self):
        airmass = np.array([5.0, 4.0, 3.0, 2.9, 2.8, 2.7, 2.6, 2.5])
        signal = 4.0 - 0.2 * airmass
        signal += np.array([0.01, -0.01, 0.015, 0.0, -0.01, 0.01, -0.005, 0.0])
        signal[4] += 1.0

        result = fit_langley(pd.Series(signal, index=airmass), robust=True, weighted=True)

        self.assertAlmostEqual(result['res_series']['slope'], -0.20334998395164727)
        self.assertAlmostEqual(result['res_series']['intercept'], 4.0151440866569414)


if __name__ == '__main__':
    unittest.main()
