from numpy.distutils.core import setup, Extension
from distutils.dir_util import copy_tree
from distutils.errors import DistutilsFileError

#import shutil


def get_version_number():
    for l in open('symgroupy/__init__.py', 'r').readlines():
        if not(l.find('__version__')):
            exec(l, globals())
            return __version__


# Make python package

try:
    copy_tree('../src', './src', update=True)
except DistutilsFileError:
    pass

s_dir = 'src/'

symgroupy = Extension('symgroupy.symgrouplib',
                      #include_dirs=include_dirs_numpy,
                      extra_f77_compile_args=['-ffixed-line-length-0'],
                      libraries=['lapack', 'blas'],
                      sources=['symgrouplib.pyf',
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
      author='Abel Carreras',
      author_email='abelcarreras83@gmail.com',
      packages=['symgroupy'],
      ext_modules=[symgroupy])
