Source code organization
============================================

*fwdpy* is a complex package made up of source files in several different languages/grammars.

The point of this section is to give some guidance to anyone who wants to modify the source.

The source code sits in two different directory:

1. fwdpy: this is where definitions are.  Files in this directory are from a mix of sources.
2. include: this is where any C++ headers containing type declarations are.

Below, all file types are fair game for modification unless otherwise noted.
   
In the fwdpy directory, you will find the following types of files:

1. .pyx/.pxd are Cython files.
2. .py files are Python files written by the developers
3. .cc files are C++ written by the developers
4. .cpp files are C++ files written by the Cython compiler.  **DO NOT MODIFY THESE**.  These files are affected by the contents of the .pyx/.pxd files

In the include directory, you fill find the following types of files:

1. .hpp are C++ headers written by the developers.

When reading the C++ source files (.cc/.hpp), headers that are part of this package are included using double quotes, meaning that you will find the file in the include directory of the source repo.
   
