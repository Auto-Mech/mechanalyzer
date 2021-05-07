""" Install Mechanism Analyzer functions
"""

from distutils.core import setup


setup(
    name="mechanalyzer",
    version="0.1.0",
    packages=[
        'mechanlyzer',
        'mechanalyzer.builder',
        'mechanalyzer.calculator',
        'mechanalyzer.parser',
        'mechanalyzer.plotter',
        'mechanalyzer.inf'
    ],
    package_dir={
          'builder': 'builder',
          'calculator': 'calculator',
          'parser': 'parser',
          'plotter': 'plotter',
          'inf': 'inf'
    },
    package_data={
        'mechanalyzer': ['tests/data/*.txt',
                         'tests/data/*.dat',
                         'tests/data/*.csv']
    }
)
