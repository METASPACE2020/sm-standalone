from setuptools import setup, find_packages
from distutils.core import Extension
from shutil import copyfile
import os

# before building the wheel, CMake must be run manually from cpp/build;
# it requires git, gcc>=4.8, boost, libxml2, and libfftw3 to be installed

import sys
sys.path.insert(0, os.path.join("python", "fastims"))
from utils import shared_lib

shared_lib_filename = shared_lib('ims_cffi')

# I couldn't make setup() to copy files in a non-trivial way, is it possible?
extra_files = { 
  os.path.join('cpp', 'build', shared_lib_filename) :
  os.path.join('python', 'fastims', shared_lib_filename),

  os.path.join('cpp', 'cffi', 'ims.h') :
  os.path.join('python', 'fastims', 'ims.h')
}

for src, dst in extra_files.items():
  copyfile(src, dst)

setup(
  name = 'fastims',
  version = '0.0.1',
  author = 'Artem Tarasov',
  author_email = 'artem.tarasov@embl.de',
  url = 'https://github.com/SpatialMetabolomics/SM_standalone',
  license = 'Apache 2.0',
  description = 'imaging mass-spec lib',
  packages = find_packages(where='python'),
  package_dir = {'':'python'},
  package_data = {'fastims': [shared_lib_filename, 'ims.h']},
  setup_requires = ['wheel>=0.27.0'],
  install_requires = ['cffi>=1.0', 'numpy'],

  # force bdist_wheel to believe it's a platform-specific wheel
  # (WHICH IT FUCKING IS, yet old hacks found on the internet
  # no longer work since 0.27.0, and the proper way is either
  # not documented or doesn't exist at all)
  ext_modules = [Extension('fastims._dummy', sources = ['dummy.c'])],

  classifiers = [
    'Development Status :: 1 - Planning',
    'Intended Audience :: Science/Research',
    'License :: OSI Approved :: Apache Software License',
    'Topic :: Scientific/Engineering :: Bio-Informatics',
  ]
)

for dst in extra_files.values():
  os.remove(dst)

