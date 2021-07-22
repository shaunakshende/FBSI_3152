#!/usr/bin/python3
if __name__ == '__main__':
  import sys
  import os
  sys.path.insert(0, os.path.abspath('config'))
  import configure
  configure_options = [
    '--download-fblaslapack',
    '--download-mpich',
    '--with-cc=gcc',
    '--with-cxx=g++',
    '--with-fc=gfortran',
    'PETSC_ARCH=arch-linux2-c-debug',
  ]
  configure.petsc_configure(configure_options)
