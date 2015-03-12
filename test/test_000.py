#!/usr/bin/env python

def test_import_root():
    '''If this fails, your ROOT environment is broken.

    Or, the user env is not setup to find the libraries from this package.

    ROOT:

    If you installed ROOT from source, by hand, source it's "thisroot.sh".

    If you are using Debian/Ubuntu broken ROOT packages set:

      export PYTHONPATH=/usr/share/python-support/root:/usr/lib/python2.7/dist-packages

    If you are using an orchestrated installation of your software
    stack, follow whatever instructions it provides.

    User env:

    Add the library installation location:

      export LD_LIBRARY_PATH=/path/to/install/lib

    '''
    import ROOT
    
