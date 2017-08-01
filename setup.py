from setuptools import setup

def readme():
    with open('README.md') as f:
        return f.read()

setup(name='genemunge',
      version='0.0',
      description='Tools for munging genomics data',
      packages=['genemunge'],
      install_requires=[
          'h5py',
          'matplotlib',
          'numpy',
          'pandas',
          'pytest',
          'seaborn',
          'tables',
          'cytoolz'
          ],
      tests_require=[
          'pytest'
      ],
      zip_safe=False)
