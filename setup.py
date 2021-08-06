import io
import re
from glob import glob
from os.path import basename
from os.path import dirname
from os.path import join
from os.path import splitext


with io.open('README.md', 'r', encoding='utf-8') as readme_file:
    readme = readme_file.read()

from setuptools import setup, find_packages

setup(
    name = 'pybedforms',
    version = '0.1.0',
    description = 'Python implementation of Bedforms 4.0 for Simulating Bedforms and Cross-Bedding',
    long_description = readme,
    long_description_content_type='text/markdown',
    url='https://github.com/andrewannex/pybedforms',
    author = "Andrew M. Annex",
    author_email = "annex@jhu.edu",
    license='BSD-3-Clause',
    packages = find_packages('src'),
    package_dir = {'': 'src'},
    py_modules=[splitext(basename(path))[0] for path in glob('src/*.py')],
    requires=['numpy', 'tqdm', 'scikit-image', 'scipy', 'matplotlib', 'mayavi', 'PyQt5'],
    keywords=['usgs','dunes','simulation','geology','sedimentology'],
    classifiers=[
        'Natural Language :: English',
        'License :: OSI Approved :: BSD License',
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
    ]

)
