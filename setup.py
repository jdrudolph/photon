from setuptools import setup, find_packages

setup(name='PHOTON',
      version='0.1',
      description='PHOsphoproteomic dissecTiOn via Networks',
      author='Jan Daniel Rudolph',
      author_email='jan.daniel.rudolph@gmail.com',
      license='MIT',
      packages=find_packages(),
      install_requires=[
          'numpy',
          'pandas',
          'scipy',
          'matplotlib',
          'statsmodels',
          'networkx',
          'flask',
	  'celery',
	  'redis',
          'goenrich',
          'pytest-runner'],
      tests_require=['pytest'])
