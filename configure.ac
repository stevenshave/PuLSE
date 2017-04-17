AC_PREREQ([2.69])
AC_INIT([PuLSE], [1.0], [stevenshave@gmail.com])
AC_CONFIG_SRCDIR([src/PuLSE.cpp])


# Checks for programs.
AC_PROG_CXX


# Checks for typedefs, structures, and compiler characteristics.
AC_CHECK_HEADER_STDBOOL

# Checks for library functions.
AC_CHECK_FUNCS([pow])

AC_CONFIG_FILES([Makefile])
AC_LANG(C++)
CXXFLAGS="-Wall -O3 -std=c++1y"
AC_HEADER_STDC
# FLag -std=c++1y
AC_MSG_CHECKING([if c++ compiler supports -std=c++1y (C++14)])
AC_COMPILE_IFELSE([AC_LANG_PROGRAM([])], [AC_MSG_RESULT([yes]); CPPUTEST_CPP14FLAG="-std=c++1y" ], [AC_MSG_RESULT([no])])

AC_OUTPUT
