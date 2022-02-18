#!python


from setuptools import setup


setup(
    name='autogromacs',
    version='0.4',
    packages=[
        'autogromacs.core',
        'autogromacs'
    ],
    platforms='any',
    license="MIT",
    install_requires=[
        "MDAnalysis"
        ],
    entry_points={
        'console_scripts': [
            "autogromacs=autogromacs.autoGromacs:main",
            "gromacsfftester=autogromacs.forceFieldCompatTester:main",
            "mdanalyze=autogromacs.mdanalysis:main",
            "plotxvg=autogromacs.parseXVG:main",
            "mdmetrics=autogromacs.simulation_metrics:main"
        ]
    },
    package_data={
        '': ["mdp/*.mdp"]
    }
)
