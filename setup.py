from distutils.core import setup, Extension
setup(name='FastNW', version='0.1',  \
      ext_modules=[Extension('FastNW', ['FastNWModule.c'])])
