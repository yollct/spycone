from setuptools import setup, find_packages
import pathlib


README = (pathlib.Path(__file__).parent / "README.md").read_text()

def setup_package():
    metadata = dict(name='spycone',
    version='0.1.3',
    description='A splicing-aware time course network enricher',
    long_description=README,
    long_description_content_type="text/markdown",
    url='https://github.com/yollct/spycone.git',
    author='Chit Tong Lio',
    author_email='ct7298@gmail.com',
    license='GPLv3',
    classifiers=[
      "Programming Language :: Python :: 3",
      "Programming Language :: Python :: 3.7",
      "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
      "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    packages=find_packages(),
    include_package_data=True,
    python_requires='>=3.7',
    install_requires=[
      'pandas>=1.0.1',
      'numpy>=1.22',
      'seaborn',
      'scikit-learn>=0.23.2',
      'scikit-learn-extra>=0.1.0',
      'networkx>=2.4',
      'matplotlib>=3.1.3',
      'gprofiler-official>=1.0.0',
      'scipy>=1.4.1',
      'statsmodels>=0.11.0',
      'tslearn>=0.5.2',
      'python-louvain',
      'plotly>=4.14.3',
      'biopython',
      'gtfparse',  
      'joblib',
      'nease'])
  
    ext_modules = []
    metadata['setup_requires']=['numpy']

    setup(**metadata)


if __name__ == "__main__":
  setup_package()
