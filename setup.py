""" Install Mechanism Analyzer functions
"""

from distutils.core import setup

setup(name="mechanalyzer",
      version="0.1.0",
      packages=['plotter',
                'calculator'],
      package_dir={
          'plotter': 'plotter',
          'calculator': 'calculator'},
      package_data={
          'plotter': ['tests/data/*.txt', 'tests/data/*.csv'],
          'chemkin': ['tests/data/*.txt', 'tests/data/*.csv']})
