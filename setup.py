""" Install Interfaces to Routines to Fit Rate Constants
"""

from distutils.core import setup

setup(name="RateFit",
      version="0.1.0",
      packages=['ratefit',
                'ratefit.calc',
                'ratefit.fit',
                'ratefit.fit.arrhenius',
                'ratefit.fit.troe'],
      package_data={'ratefit': ['tests/data/*.txt',
                                'tests/data/*.csv',
                                'fit/arrhenius/dsarrfit.mako']}
)
