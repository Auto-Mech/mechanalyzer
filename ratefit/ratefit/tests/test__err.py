""" test error
"""

# Set of rates
# Set of single rate

def test__():
    """ test ratefit.fit.err
    """

    # Calculate the sum-of-square errors and mean-average-errors
    mean_avg_err, max_avg_err = ratefit.fit.fitting_errors(
        RATEKS, fit_ks)

    assert numpy.allclose(mean_avg_err, ref_mean_avg_err)
    assert numpy.allclose(max_avg_err, ref_max_avg_err)
