#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.67])
AC_INIT([DSQSS], [1.1.17], [amasaki@issp.u-tokyo.ac.jp] ,[dsqss])
#AC_CONFIG_AUX_DIR[(config)]
AC_CONFIG_SRCDIR([src/measure_specific.cc])
AC_CONFIG_HEADERS([config.h])
AM_INIT_AUTOMAKE([no-dist-gzip dist-bzip2])

#AM_PROG_LIBTOOL
# Checks for programs.
AC_PROG_CXX
AC_PROG_CC

# Checks for libraries.

# Checks for header files.
AC_CHECK_HEADERS([stdlib.h string.h sys/time.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_INLINE

# Checks for library functions.
AC_CHECK_FUNCS([gettimeofday sqrt strchr strpbrk])

AC_CONFIG_FILES([Makefile
                 src/Makefile])
AC_OUTPUT
