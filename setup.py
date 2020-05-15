import io
import re
from glob import glob
from os.path import basename
from os.path import dirname
from os.path import join
from os.path import splitext


from setuptools import setup, find_packages

setup(
    name = 'pybedforms',
    version = '0.0.1',
    long_description = "tbd",
    url='https://github.com/andrewannex/pybedforms',
    author = "Andrew Annex",
    author_email = "annex@jhu.edu",
    license='BSD-3',
    packages = find_packages('src'),
    package_dir = {'': 'src'},
    py_modules=[splitext(basename(path))[0] for path in glob('src/*.py')],
    requires=['numpy', 'numba', 'matplotlib'],
)