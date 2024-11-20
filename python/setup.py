from setuptools import setup
from setuptools.command.build_ext import build_ext
from setuptools.command.install import install
from distutils.dir_util import copy_tree
from distutils.errors import DistutilsFileError
import sys, os
import subprocess
import pathlib

try:
    # from wheel.bdist_wheel import bdist_wheel as _bdist_wheel
    from setuptools.command.bdist_wheel import bdist_wheel as _bdist_wheel
except ModuleNotFoundError:
    _bdist_wheel = object


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

        # make compilation dir if needed
        if not os.path.exists(self.build_temp):
            os.makedirs(self.build_temp)

        print('self.build_lib:', self.build_lib)

        # define module dir to place fortran extension
        workdir = os.path.dirname(os.path.abspath(__file__))
        install_dir = pathlib.Path(workdir, 'symgroupy')

        # build with meson
        subprocess.check_call(['meson', 'setup', self.build_temp, '--prefix', install_dir])
        subprocess.check_call(['meson', 'compile', '-C', self.build_temp])
        if '--inplace' in sys.argv:
            subprocess.check_call(['meson', 'install', '-C', self.build_temp])


class InstallWithBuildExt(install):
    def run(self):

        self.build_temp = os.path.join(os.path.abspath(os.path.dirname(__file__)), 'build/temp')

        if not os.path.exists(self.build_temp):
            os.makedirs(self.build_temp)

        self.install_lib = self.install_lib.replace('purelib/', '')

        # define module dir to install fortran extension
        install_dir = pathlib.Path(self.install_lib, 'symgroupy')
        install_dir = os.path.abspath(install_dir)

        # build with meson and install
        subprocess.check_call(['meson', 'setup', self.build_temp, '--prefix', install_dir])
        subprocess.check_call(['meson', 'compile', '-C', self.build_temp])
        subprocess.check_call(['meson', 'install', '-C', self.build_temp])

        # install
        import distutils.command.install as orig
        orig.install.run(self)


class MesonBdistWheel(_bdist_wheel):
    def run(self):

        self.build_temp = os.path.join(os.path.abspath(os.path.dirname(__file__)), 'build/temp')

        if not os.path.exists(self.dist_dir):
            os.makedirs(self.dist_dir)

        # define project root dir
        workdir = os.path.dirname(os.path.abspath(__file__))

        # define distribution dir
        dist_dir = pathlib.Path(workdir, self.dist_dir)
        dist_dir = os.path.abspath(dist_dir)

        # build with meson
        subprocess.check_call(['meson', 'setup', self.build_temp, '--prefix', dist_dir])
        subprocess.check_call(['meson', 'compile', '-C', self.build_temp])

        self.root_is_pure = False
        super().run()


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
                'bdist_wheel': MesonBdistWheel,
                },
      url='https://github.com/abelcarreras/symgroup',
      classifiers=[
          "Programming Language :: Python",
          "License :: OSI Approved :: MIT License"]
      )
