import glob
from setuptools import setup, find_packages
from distutils.extension import Extension


import numpy as np


name = 'figs'
version = '0.1'


try: 
    from Cython.Distutils import build_ext
    from Cython.Build import cythonize
except ImportError:
    use_cython = False
    print("Cython not found")
else: use_cython = True

cmdclass = {}
ext_modules = []
includes = [np.get_include()]

print(includes)



if use_cython:
    print('****** Using Cython******')
    ext_modules += [
            Extension("figs.recombination",
                ["figs/recombination.pyx"],
                include_dirs=includes,
                extra_compile_args = ["-O3", 
                    "-ffast-math", 
                    "-march=native",
                    "-fopenmp" ],
                extra_link_args=['-fopenmp'], 
                ),
            ]
    cmdclass.update({'build_ext': build_ext})

metadata = {'name': name,
            'version': version,
            'ext_modules' : cythonize(ext_modules),
            'author': 'Scott Smalls',
            'packages': ['figs',
                ],
            'package_data': {'figs' : ['data/']},
            'scripts': glob.glob('scripts/*.py'),
            }


if __name__ == '__main__':
    dist = setup(**metadata)
