Installing the dependencies
**********************************

On Ununtu Linux::

  sudo apt-get -f install libgoogle-perftools-dev google-perftools libgsl0-dev libgsl0dbl

On OS X, use Homebrew_::

  brew install google-perftools gsl

On OS X, you can also install fwdpp and libsequence::

  brew install fwdpp libsequence

Please be aware that homebrew may be behind the latest fwdpp version, in which case you will need to install from source.

To install fwdpp_ and libsequence_, please follow the instructions at the home pages of those projects. (Hint: libsequence needs to be installed before fwdpp!).

You will need a *modern* C++ compiler!  This means GCC >= 4.8.2 and/or clang++ >= 3.5.

.. _Homebrew: http;//brew.sh
.. _fwdpp: http://molpopgen.github.io/fwdpp 
.. _libsequence: http://molpopgen.github.io/libsequence/
