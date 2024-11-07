project('symgrouplib', 'c', 'fortran', version: '0.1.0')

py = import('python').find_installation(pure: false)
py_dep = py.dependency()

#locate libraries LAPACK and BLAS
lapack = dependency('lapack', required: true)
blas = dependency('blas', required: true)


# fortran source files
fortran_sources = [
    'symgrouplibmodule.c',
    'src/symgrouplib.F',
    'src/radius.F',
    'src/connectivity.F',
    'src/jacobi.F',
    'src/linear_algebra.F',
    'src/mass.F',
    'src/measure.F',
    'src/operations.F',
]


# Fortran compiler option
fortran_flags = ['-ffixed-line-length-0', '-O2']

incdir_numpy = run_command(py,
  ['-c', 'import os; os.chdir(".."); import numpy; print(numpy.get_include())'],
  check : true
).stdout().strip()

incdir_f2py = run_command(py,
    ['-c', 'import os; os.chdir(".."); import numpy.f2py; print(numpy.f2py.get_include())'],
    check : true
).stdout().strip()

inc_np = include_directories(incdir_numpy, incdir_f2py)

symgrouplib = py.extension_module('symgrouplib',
  fortran_sources,
  incdir_f2py / 'fortranobject.c',
  include_directories: inc_np,
  dependencies : [py_dep, lapack, blas],
  fortran_args: fortran_flags,
  install : true,
  install_dir : get_option('prefix')
)



