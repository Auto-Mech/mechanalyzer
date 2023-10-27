""" Install Mechanism Analyzer functions
"""

from distutils.core import setup


setup(
    name="mechanalyzer",
    version="0.9.1",
    packages=[
        'mechanalyzer',
        'mechanalyzer.builder',
        'mechanalyzer.calculator',
        'mechanalyzer.parser',
        'mechanalyzer.plotter',
        'mechanalyzer.inf',
        'ratefit',
        'ratefit.calc',
        'ratefit.fit',
        'ratefit.fit.troe',
        'thermfit',
        'thermfit.cbh'
    ],
    package_dir={
        'mechanalyzer': 'mechanalyzer',
        'ratefit': 'ratefit',
        'thermfit': 'thermfit'
    },
    package_data={
        'mechanalyzer': ['tests/data/*.txt',
                         'tests/data/*.dat',
                         'tests/data/*.csv'],
        'ratefit': ['tests/data/*',
                    'fit/arrhenius/dsarrfit.mako',
                    'fit/troe/troefit.mako'],
        'thermfit': ['thermdb/*.csv']
    },
    scripts=[
        'mechanalyzer_bin/expand_species.py',
        'mechanalyzer_bin/ste_mech2.py',
        'mechanalyzer_bin/run_sort.py',
        'mechanalyzer_bin/generate_nasapoly.py',
        'mechanalyzer_bin/prompt.py',
        'mechanalyzer_bin/plot_pes_ckin.py',
        'ratefit_bin/fit_rates.py'
    ]
)
