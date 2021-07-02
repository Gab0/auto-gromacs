#!python


from setuptools import setup


setup(
    name='autogromacs',
    version='0.3',
    packages=[
        'autogromacs.core',
        'autogromacs'
    ],
    platforms='any',
    license="MIT",
    entry_points={
        'console_scripts': [
            "autogromacs=autogromacs.autoGromacs:main",
            "gromacsfftester=autogromacs.forceFieldCompatTester:main",
            "mdanalyze=autogromacs.mdanalysis:main",
            "plotxvg=autogromacs.parseXVG:main"
        ]
    },
    package_data={
        '': ["mdp/*.mdp"]
    }
)
