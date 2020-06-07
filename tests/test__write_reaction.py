"""
 tests writers
"""

import chemkin_io


def test__troe_writer():
    """ test chemkin_io.writer.reaction.troe
    """
    reaction1 = 'CH2OH=CH3O'
    high_params = [111111.1111, 1.11, 111.111]
    low_params = [222222.2222, 2.22, 222.222]
    troe_params1 = [333.333, 4.212, 123211334344.2322]
    troe_params2 = [333.333, 4.212, 123211334344.2322, 34958029.333422]
    colliders = (('H2O', 2.00), ('CO', 1.52), ('C2H2', 4.32))
    troe_str1 = chemkin_io.writer.reaction.troe(
        reaction1, high_params, low_params, troe_params1, colliders=())
    troe_str2 = chemkin_io.writer.reaction.troe(
        reaction1, high_params, low_params, troe_params1, colliders=colliders)
    troe_str3 = chemkin_io.writer.reaction.troe(
        reaction1, high_params, low_params, troe_params2, colliders=colliders)
    print('\ntroe_str1')
    print(troe_str1)
    print('\ntroe_str2')
    print(troe_str2)
    print('\ntroe_str3')
    print(troe_str3)


def test__lindemann_writer():
    """ test chemkin_io.writer.reaction.lindemann
    """
    reaction1 = 'CH2OH=CH3O'
    high_params = [111111.1111, 1.11, 111.111]
    low_params = [222222.2222, 2.22, 222.222]
    colliders = (('H2O', 2.00), ('CO', 1.52))
    lindemann_str1 = chemkin_io.writer.reaction.lindemann(
        reaction1, high_params, low_params, colliders=())
    lindemann_str2 = chemkin_io.writer.reaction.lindemann(
        reaction1, high_params, low_params, colliders=colliders)
    print('\nlindemann_str1')
    print(lindemann_str1)
    print('\nlindemann_str2')
    print(lindemann_str2)


def test__plog_writer():
    """ test chemkin_io.writer.reaction.plog
    """

    reaction1 = 'CH2+H=CH3'
    reaction2 = 'CH3+H=CH4'
    rate_params_dct1 = {
        1: [111111.1111, 1.11, 111.111],
        10: [222222.2222, 2.22, 222.222],
        100: [333333.3333, 3.33, 333.333],
        'high': [444444.4444, 4.44, 444.444]
    }
    rate_params_dct2 = {
        1: [111111.1111, 1.11, 111.111, 555555.5555, 5.55, 555.555],
        10: [222222.2222, 2.22, 222.222, 666666.6666, 6.66, 666.666],
        100: [333333.3333, 3.33, 333.333, 777777.7777, 7.77, 777.777],
        'high': [444444.4444, 4.44, 444.444, 888888.8888, 8.88, 888.888]
    }
    temp_dct = {
        1: [100, 500],
        10: [100, 700],
        100: [100, 1000],
        'high': [100, 1500]
    }
    err_dct = {
        1: [1.11, 11.111],
        10: [2.22, 22.222],
        100: [3.33, 33.333],
        'high': [4.44, 44.444]
    }

    plog_str1 = chemkin_io.writer.reaction.plog(
        reaction1, rate_params_dct1, {})
    plog_str2 = chemkin_io.writer.reaction.plog(
        reaction2, rate_params_dct2, temp_dct, err_dct)
    print('\nplog_str1')
    print(plog_str1)
    print('\nplog_str2')
    print(plog_str2)


def test__chebyshev_writer():
    """ test chemkin_io.writer.reaction.chebyshev
    """
    reaction = 'CH2OH=CH3O'
    high_params = [1.00, 0.00, 0.00]
    tmin, tmax = (300.000, 2200.000)
    pmin, pmax = (0.010, 98.702)
    alpha_matrix = [
        [8.684e+00, 7.500e-01, -7.486e-02, 1.879e-15],
        [-2.159e-01, 9.899e-02, 2.292e-02, 2.929e-17],
        [-1.557e-15, -3.331e-16, 3.324e-17, -8.346e-31],
        [2.159e-01, -9.899e-02, -2.292e-02, -2.929e-17],
        [-2.684e+00, -7.500e-01, 7.486e-02, -1.879e-15],
        [2.159e-01, -9.899e-02, -2.292e-02, -2.929e-17]
    ]

    cheb_str = chemkin_io.writer.reaction.chebyshev(
        reaction, high_params, alpha_matrix, tmin, tmax, pmin, pmax)
    print('\ncheb_str')
    print(cheb_str)


if __name__ == '__main__':
    test__troe_writer()
    test__lindemann_writer()
    test__plog_writer()
    test__chebyshev_writer()
