from numpy.distutils.core import setup, Extension
from distutils.dir_util import copy_tree
from distutils.errors import DistutilsFileError
import sys

def get_version_number():
    main_ns = {}
    for line in open('symgroupy/__init__.py', 'r').readlines():
        if not(line.find('__version__')):
            exec(line, main_ns)
            return main_ns['__version__']

# Make python package

try:
    copy_tree('../src', './src', update=True)
except DistutilsFileError:
    pass

s_dir = 'src/'

if sys.version_info[1] > 7:
    pyf_file = 'symgrouplib.pyf'
else:
    pyf_file = 'symgrouplib_legacy.pyf'


symgroupy = Extension('symgroupy.symgrouplib',
                      #include_dirs=include_dirs_numpy,
                      extra_f77_compile_args=['-ffixed-line-length-0'],
                      libraries=['lapack', 'blas'],
                      sources=[pyf_file,
                               s_dir + 'symgrouplib.F',
                               s_dir + 'radius.F',
                               s_dir + 'connectivity.F',
                               s_dir + 'jacobi.F',
                               s_dir + 'linear_algebra.F',
                               s_dir + 'mass.F',
                               s_dir + 'measure.F',
                               s_dir + 'operations.F'])

setup(name='symgroupy',
      version=get_version_number(),
      description='symgroupy',
      # long_description=open('readme.md').read(),
      # long_description_content_type='text/markdown',
      author='Abel Carreras',
      author_email='abelcarreras83@gmail.com',
      packages=['symgroupy'],
      ext_modules=[symgroupy],
      url='https://github.com/abelcarreras/symgroup',
      classifiers=[
          "Programming Language :: Python",
          "License :: OSI Approved :: MIT License"]
      )
