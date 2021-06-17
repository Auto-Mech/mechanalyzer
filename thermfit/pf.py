""" Handle arithmetic operations for partition functions and related
    information.

    only works for one value for pfs?
"""

import numpy


def combine(pfs, coeffs, operators):
    """ Combine sets of partition function data together via some set
        of input coefficients and arithmetic operations.

        _pf = [temps, logq, dq_dt, d2q_dt2]
    """

    assert len(pfs) == len(coeffs)

    for midx, _pf in enumerate(pfs):
        if midx == 0:
            coeff = coeffs[midx]
            final_pf = _pf
        else:
            pf2 = _pf
            coeff = coeffs[midx]
            operator = operators[midx-1]
            if coeff < 0:
                coeff = abs(coeff)
                if operator == 'multiply':
                    operator = 'divide'
            final_pf = _combine_pfs(final_pf, pf2, coeff, operator)

    return final_pf


def _combine_pfs(pfa, pfb, coeff, operator):
    """ Combine two sets of partition function data together under
        multiplication or division where the combined data is
        weighted by some coefficient.

        pfn = ((temp1, logq1, dqdt1, d2gdt21), (temp2, logq2, dqdt2, d2gdt22),)

        :param pfa: partition function set A
        :param pfb: partition function set B
        :param coeff: weighting coefficient
        :param operator: multiplication/division operator
    """

    tempsa, logqa, dq_dta, d2q_dt2a = pfa
    _, logqb, dq_dtb, d2q_dt2b = pfb

    if operator == 'multiply':
        logq = [a+b+numpy.log(coeff) for a, b in zip(logqa, logqb)]
        dq_dt = [a+b+numpy.log(coeff) for a, b in zip(dq_dta, dq_dtb)]
        d2q_dt2 = [a+b+numpy.log(coeff) for a, b in zip(d2q_dt2a, d2q_dt2b)]
    elif operator == 'divide':
        logq = [a-b-numpy.log(coeff) for a, b in zip(logqa, logqb)]
        dq_dt = [a-b-numpy.log(coeff) for a, b in zip(dq_dta, dq_dtb)]
        d2q_dt2 = [a-b-numpy.log(coeff) for a, b in zip(d2q_dt2a, d2q_dt2b)]

    return tempsa, tuple(logq), tuple(dq_dt), tuple(d2q_dt2)
