from numpy.distutils.core import setup, Extension


def get_version_number():
    for l in open('symgroupy/__init__.py', 'r').readlines():
        if not(l.find('__version__')):
            exec(l, globals())
            return __version__


symgroupy = Extension('symgroupy.symgrouplib',
                      #include_dirs=include_dirs_numpy,
                      extra_f77_compile_args=['-ff2c -ffixed-line-length-0'],
                      libraries=['lapack', 'blas'],
                      sources=['symgrouplib.pyf',
                               '../src/symgrouplib.F',
                               '../src/radius.F',
                               '../src/connectivity.F'])


setup(name='symgroupy',
      version=get_version_number(),
      description='symgroupy',
      author='Abel Carreras',
      author_email='abelcarreras83@gmail.com',
      packages=['symgroupy'],
      ext_modules=[symgroupy])
