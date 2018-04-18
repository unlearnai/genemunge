from setuptools import setup, find_packages
from setuptools.command.develop import develop
from setuptools.command.install import install


def setup_data():
    """
    Downloads and processes some data files used by genemunge.

    Args:
        None

    Returns:
        None

    """
    import genemunge, subprocess, os
    genemunge.data.downloads.download_everything(force=True)
    genemunge.data.parse_go.make_godict(genemunge.data.parse_go.GOFILE, force=True)
    # process the gene attributes
    genemunge.data.gene_attributes.create_attributes_file()
    # download the gtex data
    package_path = os.path.dirname(os.path.abspath(genemunge.__file__))
    script_path = os.path.join(package_path, 'data', 'gtex')
    subprocess.call(
            "Rscript " + os.path.join(script_path, 'recount.R') + ' ' + script_path,
            shell=True)
    genemunge.data.gtex.process_gtex.create_tissue_stats()
    genemunge.data.cleanup.remove_installed_data_files()


class PostDevelopCommand(develop):
    """Post-installation data download and processing for development mode."""
    def run(self):
        setup_data()
        develop.run(self)


class PostInstallCommand(install):
    """Post-installation data download and processing for installation mode."""
    def run(self):
        install.run(self)


def readme():
    with open('README.rst') as f:
        return f.read()


setup(name='genemunge',
      version='0.0',
      description='Tools for munging genomics data',
      long_description=readme(),
      url='http://github.com/unlearnai/genemunge',
      author='Unlearn.AI, Inc.',
      author_email='drckf@unlearn.ai',
      license='MIT',
      classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
      ],
      keywords='bioinformatics genomics',
      packages=find_packages(),
      package_data={'genemunge': ['data/gene_attributes.json',
                                  'data/go.json',
                                  'data/hgnc_complete_set.txt',
                                  'data/gtex/gene_info.csv',
                                  'data/gtex/tissue_stats.h5']},
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
      python_requires='>=3.5',
      zip_safe=False,
      cmdclass={
              'develop': PostDevelopCommand,
              'install': PostInstallCommand
              }
    )
