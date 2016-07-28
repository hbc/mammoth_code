from setuptools import setup, find_packages

from setuptools import setup


setup(
    name='mammoth',
    version='0.1.1',
    author='',
    author_email='',
    packages=['mammoth'],
    url='',
    install_requires = ['gffutils', 'joblib', 'toolz', 'logging', 'uuid'],
    license='See LICENSE.txt',
    description='',
     scripts=['scripts/mammoth-run.py'],
    long_description=open('README.txt').read(),
)
