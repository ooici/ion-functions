ION Functions
==============

Functions for utilization in the ION Parameter Function framework

Reference Information:  *add link*


#Prerequisites

This assumes basic development environment setup (git, directory structure). Please follow the
"New Developers Tutorial" for basic steps.


**Install the following if not yet present:**

**Install** git 1.7.7:
*Download the Mac or Linux installer and run it*

*OS Packages and package management:*
For Mac, use homebrew

    /usr/bin/ruby -e "$(curl -fsSL https://raw.github.com/gist/323731)"

  * python 2.7

**Install** python, hdf5 and netcdf with Homebrew
    
    brew install python

You can even reinstall git using brew to clean up your /usr/local directory
Be sure to read the pyon README for platform specific guidance to installing
dependent libraries and packages.
Linux: Note that many installs have much older versions installed by default.
You will need to upgrade couchdb to at least 1.1.0.

##libgswteos Dependency

A very important component is the libgswteos-10 library.  Installation is quite straightforward on Mac OSX and a little more hairy on Linux.

**On OSX**  

The libgswteos dependency is brew installable:

    brew tap lukecampbell/homebrew-libgswteos
    brew install libgswteos-10
    brew test -v libgswteos-10
    
**On Linux**

The dependencies for building/installing the library are:  **autoconf**, **automake**, & **libtool**

1. Obtain the tarball from:  https://github.com/lukecampbell/gsw-teos/tarball/v3.0r4
    * sha1: 6ae190b7da78d6aff7859e7d1a3bb027ce6cc5f3

1. Build Procedure

        bash ./autogen.sh  
        ./configure --prefix=/usr/local/libgswteos-10  
        make  
        sudo make install  

1. Linking Procedure:

        sudo ln -s /usr/local/libgswteos/lib/libgswteos-10.la /usr/local/lib/  
        sudo ln -s /usr/local/libgswteos/lib/libgswteos-10.so.3 /usr/local/lib/  
        sudo ln -s /usr/local/libgswteos/lib/libgswteos-10.so.3.0.0 /usr/local/lib/  
        sudo ln -s /usr/local/libgswteos/lib/libgswteos-10.so /usr/local/lib/  
        sudo ln -s /usr/local/libgswteos/include/gswteos-10.h /usr/local/include/  

1. Ensure that the global C_INCLUDE_PATH and LD_LIBRARY_PATH includes /usr/local/lib in all profiles otherwise python won't run correctly:

        export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib  
        export C_INCLUDE_PATH=$C_INCLUDE_PATH:/usr/local/include  

**Python packages and environment management:**

**Install** pip

    easy_install pip

**Install** virtualenv and virtualenvwrapper modules for your python 2.7 installation
*Note: This may require Mac's XCode (use XCode 3.3 free version*

    easy_install --upgrade virtualenv
    easy_install --upgrade virtualenvwrapper


Setup a virtualenv to run ion-functions (use any name you like):

    mkvirtualenv --python=python2.7 ion_functions

#Source

To obtain the ion-functions project, begin by [forking the GitHub repository](https://github.com/ooici/ion-functions/).  Next, clone your fork of the repository to your local machine (replace **your_name** with the name of your github account:

    git clone git@github.com:your_name/ion-functions.git
    cd ion-functions

#Installation
**Ensure you are in a virtualenv prior to running the steps below**

From the *ion-functions* directory, run the following commands:

    python bootstrap.py -v 1.7
    bin/buildout

Once those steps complete, you should be able to run the unit tests

#Running unit tests (via [nose](https://nose.readthedocs.org/en/latest/))

From the *coverage-model* directory, run the following command:

    bin/nosetests -v

This will run all tests in the ion-functions repository.  The **-v** flag denotes verbose output (the name of each test prints as it runs).  For more *nose* options, refer to the [nose documentation](https://nose.readthedocs.org/en/latest/man.html)

#Libraries Currently Included
* [Numpy](http://www.scipy.org/Tentative_NumPy_Tutorial) – array manipulation
    * import numpy as np
* [Numexpr](https://code.google.com/p/numexpr/) – relatively trivial "one line" expressions
    * import numexpr
        * vals = umexpr.evaluate('sin(x)**10 – y', local\_dict={'x': x\_vals, 'y': y\_vals})
* [Gibbs Seawater equations – from TEOS-10](https://pypi.python.org/pypi/pygsw) _(Contact Luke Campbell if you notice something's missing)_
    * from pygsw import vectors as gsw
        * vals = gsw.sp_from_sa(input_1, …)
