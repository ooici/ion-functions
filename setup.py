from setuptools import setup, find_packages

packages = find_packages()

classifiers = ''' Intended Audience :: Science/Research
Intended Audience :: Developers
Intended Audience :: Education
Operating System :: OS Independent
Programming Language :: Python
Topic :: Scientific/Engineering
Topic :: Education
Topic :: Software Development :: Libraries :: Python Modules'''

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
        setup_requires=[],
        )




