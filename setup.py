
from setuptools import setup, find_packages
from distutils.extension import Extension

import numpy as np
# setuptools DWIM monkey-patch madness
# http://mail.python.org/pipermail/distutils-sig/2007-September/thread.html#8204
import sys

if 'setuptools.extension' in sys.modules:
    m = sys.modules['setuptools.extension']
    m.Extension.__dict__ = m._Extension.__dict__


packages = find_packages()

classifiers = ''' Intended Audience :: Science/Research
Intended Audience :: Developers
Intended Audience :: Education
Operating System :: OS Independent
Programming Language :: Python
Topic :: Scientific/Engineering
Topic :: Education
Topic :: Software Development :: Libraries :: Python Modules'''


qc_extension = Extension("ion_functions.qc.qc_extensions", ["ion_functions/qc/qc_extensions.pyx", "extensions/stuck.c"], include_dirs=[np.get_include(), "extensions/"])


setup(name = 'ion-functions', 
        version='0.0.1',
        description='Python Function collection for ION',
        long_description=open('README.md').read(),
        license='LICENSE.txt',
        author='Luke Campbell',
        author_email='lcampbell@asascience.com',
        url='https://github.com/ooici/ion-functions/',
        classifiers=classifiers.split('\n'),
        packages=packages,
        keywords=['oceanography', 'seawater'],
        ext_modules=[qc_extension],
        setup_requires=['setuptools_cython'],
        install_requires=[
            'ipython==0.13.0',
            'readline',
            'numexpr==2.1',
            'nose==1.1.2',
            'pygsw==0.0.9',
            'geomag==0.9',
            'scipy==0.11.0',
            'cython'
        ]
)

