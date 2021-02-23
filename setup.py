#!python


from setuptools import setup


setup(
    name='autogromacs',
    version='0.1',
    packages=[
        'autogromacs.core',
        'autogromacs'
    ],
    platforms='any',
    license="MIT",
    entry_points={
        'console_scripts': [
            "autogromacs=autogromacs.autoGromacs:main"
        ]
    },
	package_data={
	'': ["mdp/*.mdp"]}


)
