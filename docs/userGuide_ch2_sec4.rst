License, Distribution and Installation
======================================

simuPOP is distributed under a GPL license and is hosted at\
``http://simupop.sourceforge.net``, the world's largest development and download
repository of Open Source code and applications. simuPOP is available on any
platform where Python is available, and is currently tested under both 32 and 64
bit versions of Windows (Windows 2000 and later), Linux (Redhat and Ubuntu),
MacOS X and Sun Solaris systems. Different C++ compilers such as Microsoft
Visual C++, gcc and Intel icc are supported under different operating systems.
Standard installation packages are provided for Windows, Linux, and MacOS X
systems.

If a binary distribution is unavailable for a specific platform, it is usually
easy to compile simuPOP from source, following the standard ````python setup.py
install''`` procedure. Please refer to the ``installation`` section of the
simupop website for instructions for specific platforms and compilers.

simuPOP is available for Python 2.4 and later, including the new Python 3.x
releases. Although Python 3 is incompatible with Python 2 in many ways, examples
in this guide are written in a style that is compatible with both versions of
Python. Some non-classic usages include the use of ``a//b`` instead of ``a/b``
for floored division and ``list(range(3))`` instead of ``range(3)`` for sequece
``[0,1,2]`` In particular, we use   ::

   print("Population size is %d" % size)

instead of   ::

   print "Population size is %d" % size

to output strings because the former is valid in Python 2.x (print a tuple with
one element) and will generate the same output in Python 3.x. Of course, users
of simuPOP can choose to use other styles.

Thanks to the 'glue language' nature of Python, it is easy to inter-operate with
other applications within a simuPOP script. For example, users can call any R
function from Python/simuPOP for the purposes of visualization and statistical
analysis, using ``R`` and a Python module ``RPy``. Because simuPOP utility
modules such as ``simuPOP.plotter`` and :mod:`simuPOP.sampling` makes use of
``R`` and ``rpy`` (not ``rpy2``) to plot figures, **it is hihgly recommended
that you install R and RPy with simuPOP**. In addition, although simuPOP uses
the standard ``Tkinter`` GUI toolkit when a graphical user interface is needed,
it can make use of a ``wxPython`` toolkit if it is available.


