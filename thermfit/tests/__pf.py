""" test thermfit.pf
"""

import thermfit.pf


# [temps, logq, dq_dt, d2q_dt2]
PFS = (
    # Partition Function Set 1
    ((500.0, 67.50302, 0.040505, 4.23039e-05),
     (1000.0, 75.20304, 0.0354286, 5.36839e-04),
     (1500.0, 84.2003, 0.0673282, 7.793224e-04),
     (2000.0, 103.40404, 0.0937334, 9.837278e-04)),
    # Partition Function Set 2
    ((500.0, 65.9811, 0.010705, 1.21039e-06),
     (1000.0, 72.3157, 0.0158886, 1.70879e-05),
     (1500.0, 82.8571, 0.0271782, 2.79324e-05),
     (2000.0, 100.375, 0.0437624, 3.83778e-05))
)

COEFFS = (0.5, 1.2)
OPERATORS = ('multiply', 'divide')


def test__combine():
    """ test thermfit.pf.combine
    """

    combined_pfs = thermfit.pf.combine(PFS, COEFFS, OPERATORS)
    print(combined_pfs)


if __name__ == '__main__':
    test__combine()
