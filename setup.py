# -*- coding: utf-8 -*-

from setuptools import setup, find_packages


with open('README.rst') as f:
    readme = f.read()

with open('LICENSE') as f:
    license = f.read()

setup(
    name='falcon-tools',
    version='0.0.1',
    description='Tool to inspect consensus errors in falcon assemblies',
    long_description=readme,
    author='Greg Concepcion',
    author_email='gconcepcion@pacificbiosciences.com',
    url='https://github.com/gconcepcion/falcon-probe',
    license=license,
    packages=find_packages(exclude=('tests', 'docs')),
    scripts=['falcon_tools/plot_distributions.py']
    #entry_points={'console_scripts': ['falcon_probe = falcon_probe.cli:main']}
)

