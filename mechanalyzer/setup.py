""" Install Mechanism Analyzer functions
"""

from distutils.core import setup


setup(name="mechanalyzer",
      version="0.1.0",
      packages=[
            'plotter',
            'builder',
            'calculator',
            'parser',
            'plotter',
            'par',
            'inf'],
      package_dir={
          'plotter': 'plotter',
          'builder': 'builder',
          'calculator': 'calculator',
          'parser': 'parser',
          'par': 'par',
          'inf': 'inf'},
      package_data={
          'plotter': ['tests/data/*.txt', 'tests/data/*.csv'],
          'chemkin': ['tests/data/*.txt', 'tests/data/*.csv']})
