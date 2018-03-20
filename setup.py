from setuptools import setup
from setuptools.command.develop import develop
from setuptools.command.install import install


class PostDevelopCommand(develop):
    """Post-installation data download and processing for development mode."""
    def run(self):
        import genemunge
        genemunge.data.downloads.download_everything(force=True)
        genemunge.data.parse_go.make_godict(genemunge.data.parse_go.GOFILE, force=True)
        develop.run(self)


class PostInstallCommand(install):
    """Post-installation data download and processing for installation mode."""
    def run(self):
        import genemunge
        genemunge.data.downloads.download_everything(force=True)
        genemunge.data.parse_go.make_godict(genemunge.data.parse_go.GOFILE, force=True)
        install.run(self)


def readme():
    with open('README.md') as f:
        return f.read()


setup(name='genemunge',
      version='0.0',
      description='Tools for munging genomics data',
      long_description=readme(),
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
      zip_safe=False,
      cmdclass={
              'develop': PostDevelopCommand,
              'install': PostInstallCommand
              }
    )
