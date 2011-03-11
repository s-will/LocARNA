#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT([LocARNA], [1.6.1], [will@informatik.uni-freiburg.de], [locarna],
        [http://www.bioinf.uni-freiburg.de/Software/LocARNA/])

# Define here to keep all version information in one place.
AC_SUBST([SO_VERSION], [1:0:0])
AC_SUBST([API_VERSION], [1.0])

AC_CONFIG_MACRO_DIR([m4])

LT_INIT

AM_INIT_AUTOMAKE

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC

AC_LANG_CPLUSPLUS

AC_PATH_PROGS([PERL],[perl perl5.004 perl5.003 perl5.002 perl5.001 perl5])

AC_PATH_PROG([HELP2MAN],[help2man],false)
AC_PATH_PROG([POD2MAN],[pod2man],false)




# Vienna RNA package library path support, if not installed in standard directory/path
AC_ARG_WITH([vrna],
    AC_HELP_STRING(
        [--with-vrna=PREFIX],
        [Installation directory of the Vienna RNA library]
      )
) 

AC_ARG_VAR([RNAfold],[RNAfold name (def=RNAfold)])
AC_ARG_VAR([RNAplfold],[RNAplfold name (def=RNAplfold)])
AC_ARG_VAR([RNAalifold],[RNAalifold name (def=RNAalifold)])


if test "$with_vrna" = "" ; then
   AC_PATH_PROG([RNAfold],[RNAfold],[false])
   AC_PATH_PROG([RNAplfold],[RNAplfold],[false])
   AC_PATH_PROG([RNAalifold],[RNAalifold],[false])
else
   AC_PATH_PROG([RNAfold],[RNAfold],[false],[$with_vrna/bin])
   AC_PATH_PROG([RNAplfold],[RNAplfold],[false],[$with_vrna/bin])
   AC_PATH_PROG([RNAalifold],[RNAalifold],[false],[$with_vrna/bin])  
fi


# error output if ViennaRNA not found
VRNA_OK=true;
if test "$RNAfold" = "false" -o "$RNAplfold" = "false" -o "$RNAalifold" = "false" ; then
    VRNA_OK=false
fi



# Checks for libraries.

# Checks for header files.
AC_HEADER_STDC

# which T-Coffee to use
AC_MSG_CHECKING([where to find t-coffee])
if test "$TCOFFEE" = "" ; then
  TCOFFEE=t_coffee
fi
AC_ARG_VAR(TCOFFEE,[tcoffee location or name (def=t_coffee)])
AC_MSG_RESULT([$TCOFFEE])

AC_PATH_PROG([TCOFFEE],[$TCOFFEE],false)



##################
# Debug option
#
AC_MSG_CHECKING([whether to build with debug information])
debuger=no
AC_ARG_ENABLE([debug],
    AC_HELP_STRING(
        [--enable-debug],
        [enable debug data generation (def=no)]
    ),
    debuger="$enableval"
)
AC_MSG_RESULT($debuger)
if test x"$debuger" = x"yes"; then
    AC_DEFINE([DEBUG])
    CPPFLAGS="$CPPFLAGS -g -Wall"
else
    AC_DEFINE([NDEBUG])
fi

CPPFLAGS="$CPPFLAGS -Wno-deprecated"


####################
# Static linking
#
AC_MSG_CHECKING([whether to link statically])
static_linkage=no
AC_ARG_ENABLE([static],
    AC_HELP_STRING(
        [--enable-static],
        [link statically (def=no)]
    ),
    static_linkage="$enableval"
)
AC_MSG_RESULT($static_linkage)
if test x"$static_linkage" = x"yes"; then
   LDFLAGS="$LDFLAGS -static"
fi


## flag for use of long double partition function type
AC_MSG_CHECKING([whether to support large partition function])
use_large_pf=no
AC_ARG_ENABLE([large-pf],
    AC_HELP_STRING(
        [--enable-large-pf],
        [support very large partition functions (def=no)]
    ),
    use_large_pf="$enableval"
)
AC_MSG_RESULT($use_large_pf)
if test x"$use_large_pf" = x"yes"; then
   AC_DEFINE([LARGE_PF])
fi



# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_CONST
AC_TYPE_SIZE_T

AC_CONFIG_FILES([Makefile LocARNA-${API_VERSION}.pc:LocARNA.pc.in])


# AC_FUNC_MALLOC


AC_OUTPUT([src/Makefile Utils/Makefile lib/Makefile lib/perl/Makefile lib/perl/MLocarna/Makefile Utils/mlocarna Utils/locarnate Utils/aln2fa.pl])


if test "$VRNA_OK" = "false" ; then
   AC_MSG_NOTICE([============================================================])
   AC_MSG_NOTICE()
     if test "$with_vrna" != "" ; then
        AC_MSG_NOTICE(Cannot find Vienna RNA package in given path '$with_vrna'.)
     else
        AC_MSG_NOTICE(Cannot find Vienna RNA package in standard path, for non-standard path use --with-vrna=PREFIX.)
     fi
    AC_MSG_NOTICE()
    AC_MSG_NOTICE(The Vienna RNA package is required for full functionality (use version >= 1.8.*).)
    AC_MSG_NOTICE(It can be obtained from http://www.tbi.univie.ac.at/~ivo/RNA/.)
    AC_MSG_NOTICE()
fi
