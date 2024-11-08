from setuptools import setup
from setuptools.command.build_ext import build_ext
from setuptools.command.install import install
from distutils.dir_util import copy_tree
from distutils.errors import DistutilsFileError
import sys, os
import subprocess
import pathlib

try:
    from wheel.bdist_wheel import bdist_wheel as _bdist_wheel
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


from setuptools.command.egg_info import egg_info as _egg_info

class CustomEggInfo(_egg_info):
    """ Custom egg_info command to set the correct egg-info name and directory """

    def run(self):
        # Defineix el directori on es guardarà l'egg-info
        egg_info_dir = os.path.join(os.getcwd(), 'build', 'egg-info')
        egg_info_dir = '/Users/abel/Programes/symgroup/python/build/'
        print('anterior', self.egg_base)

        # self.distribution.package_dir = egg_info_dir

        self.egg_base = egg_info_dir
        self.egg_base = '/Users/abel/Programes/symgroup/python/build/bdist.macosx-10.9-universal2/'
        print(f"Guardant egg-info a {egg_info_dir}")
        #self.do_egg_install()

        print('egg_base: ',self.egg_base)
        print('egg_name: ',self.egg_name)
        print('egg_info: ',self.egg_info)
        print('egg_version: ',self.egg_version)
        print('broken_egg_info: ', self.broken_egg_info)

        self.egg_info = self.egg_base + 'egg/' + self.egg_name + '-' + self.egg_version + '-' + 'py3.11' + '.egg-info'
        print('self.egg_info: ', self.egg_info)

        self.egg_info = 'build/bdist.macosx-10.9-universal2/egg/symgroupy-0.5.12-py3.11.egg-info'

        # ValueError: symgroupy-0.5.12-py3.11.egg-info but not found
        # Crida a la funció base per generar l'egg-info correctament format
        super().run()


class MesonBuildExt(build_ext):
    def run(self):
        # make compilation dir
        if not os.path.exists(self.build_temp):
            os.makedirs(self.build_temp)

        print('self.build_lib:', self.build_lib)

        workdir = os.path.dirname(os.path.abspath(__file__))

        install_dir = pathlib.Path(workdir, 'symgroupy')
        # install_dir = os.path.abspath(install_dir)

        #subprocess.check_call([sys.executable, 'setup.py', 'egg_info'])

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
        print('Install dir: ', self.install_lib, install_dir)
        print('build dir:', self.build_temp)
        install_dir = os.path.abspath(install_dir)

        print('version', self)
        # subprocess.check_call([sys.executable, 'setup.py', 'egg_info', '--egg-base', self.install_lib])
        #subprocess.check_call([sys.executable, 'setup.py', 'egg_info'])

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

        workdir = os.path.dirname(os.path.abspath(__file__))

        dist_dir = pathlib.Path(workdir, self.dist_dir)
        dist_dir = os.path.abspath(dist_dir)

        print('dist_dir: ', dist_dir)

        subprocess.check_call(['meson', 'setup', self.build_temp, '--prefix', dist_dir])
        subprocess.check_call(['meson', 'compile', '-C', self.build_temp])
        # subprocess.check_call(['meson', 'install', '-C', self.build_temp])

        # impl_tag, abi_tag, plat_tag = self.get_tag()
        # archive_basename = "{}-{}-{}-{}".format(self.wheel_dist_name, impl_tag, abi_tag, plat_tag)
        # print('archive_basename: ', archive_basename)

        # self.skip_build = False
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
                # 'egg_info': CustomEggInfo
                },
      url='https://github.com/abelcarreras/symgroup',
      classifiers=[
          "Programming Language :: Python",
          "License :: OSI Approved :: MIT License"]
      )
