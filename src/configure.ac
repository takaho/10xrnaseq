#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([10xtools], [0.1], [takaho.endo@riken.jp])
AM_INIT_AUTOMAKE([foreign])
AC_CONFIG_SRCDIR([count_and_rename.cxx])
AC_CONFIG_HEADERS([config.h])

# Checks for programs.
AC_PROG_CXX

# Checks for libraries.
# FIXME: Replace `main' with a function in `-lboost_iostreams':
AC_CHECK_LIB([boost_iostreams], [main])
# FIXME: Replace `main' with a function in `-lboost_program_options':
AC_CHECK_LIB([boost_program_options], [main])
# FIXME: Replace `main' with a function in `-lboost_regex':
AC_CHECK_LIB([boost_regex], [main])

# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.
AC_CHECK_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_SIZE_T
AC_TYPE_UINT64_T

# Checks for library functions.

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
