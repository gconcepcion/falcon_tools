# -*- coding: utf-8 -*-

from setuptools import setup, find_packages


with open('README.rst') as f:
    readme = f.read()

with open('LICENSE') as f:
    license = f.read()

setup(
    name='falcon-tools',
    version='0.0.1',
    description='Some miscellaneous falcon tool scripts written in python',
    long_description=readme,
    author='Greg Concepcion',
    author_email='gconcepcion@pacificbiosciences.com',
    url='https://github.com/gconcepcion/falcon-probe',
    license=license,
    packages=find_packages(exclude=('tests', 'docs')),
    install_requires=[
                  'pandas==0.20.3',
                  'matplotlib==2.0.2',
                  'pbcore'
                        ],
    scripts=['bin/plot_distributions.py', 'bin/clean_fasta.py', 'bin/get_homologs.py']
    #entry_points={'console_scripts': ['falcon_probe = falcon_probe.cli:main']}
)

