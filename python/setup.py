from numpy.distutils.core import setup, Extension


symgroupy = Extension('symgroupy.symgrouplib',
                      #extra_compile_args=['-ffixed-line-length-0'],
                      #include_dirs=include_dirs_numpy,
                      #include_dirs=['../include'],
                      extra_f77_compile_args=['-ff2c -ffixed-line-length-0'],
                      libraries=['lapack', 'blas'],
                      sources=['symgrouplib.pyf',
                               '../src/symgrouplib.F',
                               '../src/radius.F',
                               '../src/connectivity.F'])


setup(name='symgroupy',
      version='0.1',
      description='symgroupy',
      author='Abel Carreras',
      author_email='abelcarreras83@gmail.com',
      packages=['symgroupy'],
      ext_modules=[symgroupy])
