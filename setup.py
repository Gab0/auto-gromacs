#!python


from setuptools import setup


setup(
    name='autogromacs',
    version='0.2',
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
            "mdanalyze=autogromacs.mdanalysis:main"
        ]
    },
	package_data={
	'': ["mdp/*.mdp"]}


)
