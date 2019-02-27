Description
------------

The toolbox contains several Matlab scripts, MEX source code in c++
and precompiled mex-files for Linux and Windows. See "contents.m" for
short description of the Matlab functions. A short demo is provided in
"demopcamv.m".

Matlabs scripts is enough for working with data matrices in which NaNs
represent missing values requires no installation. Mex-files are
needed for dealing with sparse matrices containing only observed
values.


Installation
------------

Working with data matrices in which NaNs represent missing values
requires no installation: Just add the toolbox directory into the
Matlab path.

Working with sparse matrices containing only observed values requires
building several mex-files. This distribution contains precompiled
mex-files for Linux and Windows which should work.

In case you want to build the mex-files yourselves: Typing "install"
in Matlab from the toolbox directory will create mex files with no
support of multiple threads. Typing "install_threads" creates mex
files which support multiple threads using the pthreads library.

Building mex files in Windows requires an installed c++ compiler.
Using MinGW (www.mingw.org) is a good option.

Building mex-files which support multiple threads requires the Posix
Threads library installed. In Windows, Posix Threads for Win32
(sourceware.org/pthreads-win32/) should work nicely. Modify
install_threads.m such that the path to pthreads.h is included in the
#include search and link a relevant lib file from the pthreads
package. When the functions are called in Matlab, make sure that a
necessary dll file for the pthreads library is in one of the
directories in the Windows path (check PATH environment variable).

This distribution includes Linux mex-files which support threads.  The
default Windows dlls included in this package do not support threads.
The dll-files in the "pthreads-win32" directory are build for Windows
using "libpthreadGCE2.a" and "pthreadGCE2.dll" from the Posix Threads
package for Win32.

----------------------------------------------------------------------

This software is provided "as is", without warranty of any kind.
Alexander Ilin, Tapani Raiko