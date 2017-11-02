from setuptools import setup, find_packages

import os
def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(name='photon_ptm',
        version='0.2.0',
        description='PHOsphoproteomic dissecTiOn via Networks',
        long_description=read('README.md'),
        url='http://www.github.com/jdrudolph/photon',
        author='Jan Rudolph',
        author_email='jan.daniel.rudolph@gmail.com',
        license='MIT',
        packages=find_packages(),
        install_requires=[
            'numpy',
            'pandas>=0.21',
            'scipy',
            'matplotlib',
            'statsmodels',
            'networkx',
            'perseuspy==0.3.0',
            'flask',
            'celery',
            'redis',
            'goenrich',
            'pytest-runner'],
        tests_require=['pytest'])
