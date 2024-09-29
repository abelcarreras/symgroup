from setuptools import setup
from setuptools.command.build_ext import build_ext
from setuptools.command.install import install
from distutils.dir_util import copy_tree
from distutils.errors import DistutilsFileError
import sys, os
import subprocess
import pathlib

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


class MesonBuildExt(build_ext):
    def run(self):
        # make compilation dir
        if not os.path.exists(self.build_temp):
            os.makedirs(self.build_temp)

        print('self.build_lib:', self.build_lib)

        workdir = os.path.dirname(os.path.abspath(__file__))

        install_dir = pathlib.Path(workdir, 'symgroupy')

        subprocess.check_call(['meson', 'setup', self.build_temp, '--prefix', install_dir])
        subprocess.check_call(['meson', 'compile', '-C', self.build_temp])
        if '--inplace' in sys.argv:
            subprocess.check_call(['meson', 'install', '-C', self.build_temp])


class InstallWithBuildExt(install):
    def run(self):

        self.build_temp = os.path.join(os.path.abspath(os.path.dirname(__file__)), 'build/temp')

        if not os.path.exists(self.build_temp):
            os.makedirs(self.build_temp)

        install_dir = pathlib.Path(self.install_lib, 'symgroupy')

        subprocess.check_call(['meson', 'setup', self.build_temp, '--prefix', install_dir])
        subprocess.check_call(['meson', 'compile', '-C', self.build_temp])
        subprocess.check_call(['meson', 'install', '-C', self.build_temp])


"""
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
"""

setup(name='symgroupy',
      version=get_version_number(),
      description='symgroupy',
      # long_description=open('readme.md').read(),
      # long_description_content_type='text/markdown',
      author='Abel Carreras',
      author_email='abelcarreras83@gmail.com',
      packages=['symgroupy'],
      ext_modules=[],
      cmdclass={'build_ext': MesonBuildExt,
                'install': InstallWithBuildExt,
                },
      url='https://github.com/abelcarreras/symgroup',
      classifiers=[
          "Programming Language :: Python",
          "License :: OSI Approved :: MIT License"]
      )
