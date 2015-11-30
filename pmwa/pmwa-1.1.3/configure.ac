#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.63])
AC_INIT([pmwa], [1.1.1], [akiko.masaki@riken.jp])
AC_CONFIG_SRCDIR([src/lattgene.cc])
AC_CONFIG_HEADERS([config.h])
AM_INIT_AUTOMAKE
# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AC_PROG_INSTALL
AM_PROG_CC_C_O
AC_PROG_RANLIB
# Checks for libraries.

# Checks for header files.
AC_CHECK_HEADERS([float.h stdlib.h string.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_INLINE

# Checks for library functions.
AC_CHECK_FUNCS([sqrt])

AC_CONFIG_FILES([Makefile
                 lib/Makefile
                 src/Makefile])
AC_OUTPUT
