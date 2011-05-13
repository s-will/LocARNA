#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT([LocARNA], [1.6.2pre], [will@informatik.uni-freiburg.de], [locarna],
        [http://www.bioinf.uni-freiburg.de/Software/LocARNA/])

# Define here to keep all version information in one place.
AC_SUBST([SO_VERSION], [1:0:1])
AC_SUBST([API_VERSION], [1.0])

AC_CONFIG_MACRO_DIR([m4])

LT_INIT

AM_INIT_AUTOMAKE

# Checks for programs.
AC_PROG_CXX
AC_PROG_CXXCPP
AC_PROG_CC
AC_PROG_CPP

AC_LANG_CPLUSPLUS

AC_PATH_PROGS([PERL],[perl perl5.004 perl5.003 perl5.002 perl5.001 perl5])

AC_PATH_PROG([HELP2MAN],[help2man],false)
AC_PATH_PROG([POD2MAN],[pod2man],false)



# Vienna RNA package library path support, if not installed in standard directory/path
VRNA_OK=true;

AC_ARG_WITH([vrna],
    AC_HELP_STRING(
        [--with-vrna=PREFIX],
        [Installation directory of the Vienna RNA library]
      )
)

####################
## whether to link to RNAlib and offer corresponding features
AC_MSG_CHECKING([whether to link Vienna RNA library])
ENABLE_LIBRNA=no
AC_ARG_ENABLE([librna],
    AC_HELP_STRING(
        [--enable-librna],
        [enable linking to librna and corresponding features (def=disabled)]
    ),
    ENABLE_LIBRNA="$enableval"
)
AC_MSG_RESULT([$ENABLE_LIBRNA])


# register Vienna programs as configuration variables
AC_ARG_VAR([RNAfold],[RNAfold name (def=RNAfold)])
AC_ARG_VAR([RNAplfold],[RNAplfold name (def=RNAplfold)])
AC_ARG_VAR([RNAalifold],[RNAalifold name (def=RNAalifold)])

# check for required Vienna programs
if test "$with_vrna" = "" ; then
   AC_PATH_PROG([RNAfold],[RNAfold],[false])
   AC_PATH_PROG([RNAplfold],[RNAplfold],[false])
   AC_PATH_PROG([RNAalifold],[RNAalifold],[false])
else
   AC_PATH_PROG([RNAfold],[RNAfold],[false],[$with_vrna/bin])
   AC_PATH_PROG([RNAplfold],[RNAplfold],[false],[$with_vrna/bin])
   AC_PATH_PROG([RNAalifold],[RNAalifold],[false],[$with_vrna/bin])  
fi

# trigger error if ViennaRNA programs not found
if test "$RNAfold" = "false" -o "$RNAplfold" = "false" -o "$RNAalifold" = "false" ; then
    VRNA_OK=false
fi

# if we want linking to librna
if test "$ENABLE_LIBRNA" = "yes" ; then
   # if --with-vrna given, add to flags for CPP and LD
   if test "$with_vrna" != "" ; then
      LDFLAGS="-L$with_vrna/lib $LDFLAGS"
      CPPFLAGS="-I$with_vrna/include $CPPFLAGS"
   fi
   
   # Check for Vienna RNA library (primitive check for one header) 
   AC_CHECK_HEADER("ViennaRNA/fold.h",
		   [],
 		   [VRNA_OK=false])

   ## link the Vienna library
   LIBS="-lRNA $LIBS"
   
   AC_DEFINE([HAVE_LIBRNA],[],[Whether librna is linked and available.])
fi


# ------------------------------
# handle open mp support
#  this is necessary for linking against the vienna lib 2.0
AC_OPENMP
AC_SUBST(OPENMP_CXXFLAGS)

CPPFLAGS="$OPENMP_CXXFLAGS $CPPFLAGS"
LDFLAGS="$OPENMP_CXXFLAGS $LDFLAGS"


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


### ----------------------------------------
### Test if Vienna lib available and ok
## don't allow compiling otherwise 
#
if test "$VRNA_OK" = "false" ; then
   AC_MSG_NOTICE([=================================================================])
   AC_MSG_NOTICE([ERROR:])
   AC_MSG_NOTICE()
   
   if test "$with_vrna" != "" ; then
      AC_MSG_NOTICE([Cannot find Vienna RNA package in given path '$with_vrna'.])
   else
      AC_MSG_NOTICE([Cannot find Vienna RNA package in standard path. For using a])
      AC_MSG_NOTICE([non-standard path, please specify with --with-vrna=PREFIX.])
   fi
   
   AC_MSG_NOTICE()
   AC_MSG_NOTICE([The Vienna RNA package (version >= 1.8) is required for compiling.])
   AC_MSG_NOTICE([It is freely available at http://www.tbi.univie.ac.at/~ivo/RNA/.])
   AC_MSG_NOTICE()
    
   exit -1;
fi


AC_OUTPUT([src/Makefile Utils/Makefile man/Makefile lib/Makefile lib/perl/Makefile lib/perl/MLocarna/Makefile Utils/mlocarna Utils/locarnate Utils/aln2fa.pl Utils/pp2dot])
