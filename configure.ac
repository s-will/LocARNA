#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT([LocARNA], [1.7.6], [will@informatik.uni-freiburg.de], [locarna],
        [http://www.bioinf.uni-freiburg.de/Software/LocARNA/])

# Define here to keep all version information in one place.
AC_SUBST([SO_VERSION], [1:0:1])
AC_SUBST([API_VERSION], [1.1])

AC_CONFIG_MACRO_DIR([m4])

LT_INIT                     # for libtool
AC_PROG_LIBTOOL             # deprecated; left for backward compatibility
# AC_SUBST(LIBTOOL_DEPS)

AM_INIT_AUTOMAKE

# Checks for programs.
AC_PROG_CXX
AC_PROG_CXXCPP
AC_PROG_CC
AC_PROG_CPP

AC_LANG_CPLUSPLUS

AC_PATH_PROGS([PERL],[perl perl5])

AC_PATH_PROG([HELP2MAN],[help2man],false)
AC_PATH_PROG([POD2MAN],[pod2man],false)


# Vienna RNA package library path support, if not installed in standard directory/path
VRNA_OK=true
VRNALIB_OK=true

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

AC_SUBST(LDFLAGS,[""])
AC_SUBST(CPPFLAGS,[""])
AC_SUBST(LIBS,[""])
AC_SUBST(LIBDEFS,[""]) # defines that have to be set when using the LocARNA library


# set PKG_CONFIG_PATH
if test "$with_vrna" != "" ; then
   export PKG_CONFIG_PATH="$with_vrna/lib/pkgconfig:$PKG_CONFIG_PATH"
fi

AS_IF([test "$ENABLE_LIBRNA" = "yes"],[
    PKG_CHECK_MODULES([VRNA],[RNAlib2 >= 2.0.7],[],[VRNALIB_OK=false])
    AC_DEFINE([HAVE_LIBRNA])
    LDFLAGS="$VRNA_LIBS $LDFLAGS"

    if test "$with_vrna" != "" ; then
        VRNA_CFLAGS="$VRNA_CFLAGS -I$with_vrna/include"
    fi
    CPPFLAGS="$VRNA_CFLAGS $CPPFLAGS"
])
AM_CONDITIONAL([HAVE_LIBRNA],[test "$ENABLE_LIBRNA" = "yes"])


####################
## whether to build locarna_X / ExpARNA-P functionality

AC_MSG_CHECKING([whether to build components for ExpARNA-P and SPARSE functionality])
AC_ARG_WITH([sparse],
    AC_HELP_STRING(
        [--with-sparse],
        [Build components relying on extra sparsification.]
    ),
    [with_sparse="yes"],
    [with_sparse="no"]
)
AC_MSG_RESULT([$with_sparse])

if test "$with_sparse $ENABLE_LIBRNA" = "yes no" ; then
   AC_MSG_NOTICE([=================================================================])
   AC_MSG_NOTICE([ERROR: Building 'with-sparse' requires linking to VRNA (enable-librna)])
   AC_MSG_NOTICE()
   exit -1;
fi

AM_CONDITIONAL([BUILD_SPARSE],[test "$with_sparse" = "yes"])



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
debug_it=no
AC_ARG_ENABLE([debug],
    AC_HELP_STRING(
        [--enable-debug],
        [enable debug data generation (def=no)]
    ),
    debug_it="$enableval"
)
AC_MSG_RESULT([$debug_it])
if test x"$debug_it" = x"yes"; then
    AC_DEFINE([DEBUG])
    CPPFLAGS="$CPPFLAGS -g -Wall"
else
    AC_DEFINE([NDEBUG])
fi

## CPPFLAGS="$CPPFLAGS -Wno-deprecated"



####################
# Static linking
#
AC_MSG_CHECKING([whether to link statically])
link_statically=no
AC_ARG_ENABLE([link-static],
    AC_HELP_STRING(
        [--enable-link-static],
        [link statically (def=no)]
    ),
    link_statically="$enableval"
)

AC_MSG_RESULT([$link_statically])
if test x"$link_statically" = x"yes"; then
   LDFLAGS="$LDFLAGS -all-static"
fi


## flag for use of long double partition function type
## attention: name clash with Vienna package macro LARGE_PF ==> define VERY_LARGE_PF
AC_MSG_CHECKING([whether to support large partition function])
use_large_pf=no
AC_ARG_ENABLE([large-pf],
    AC_HELP_STRING(
        [--enable-large-pf],
        [support very large partition functions (def=no)]
    ),
    use_large_pf="$enableval"
)
AC_MSG_RESULT([$use_large_pf])
if test x"$use_large_pf" = x"yes"; then
    LIBDEFS="$LIBDEFS -DVERY_LARGE_PF"
    AC_DEFINE([VERY_LARGE_PF])
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
fi

if test "$VRNALIB_OK" = "false" ; then
   AC_MSG_NOTICE([=================================================================])
   AC_MSG_NOTICE([ERROR:])
   AC_MSG_NOTICE()
   
   if test "$with_vrna" != "" ; then
      AC_MSG_NOTICE([Cannot find Vienna RNA package pkgconfig file RNAlib2.pc in '$with_vrna/lib/pkgconfig'.])
   else
      AC_MSG_NOTICE([Cannot find Vienna RNA package pkgconfig file RNAlib2.pc in default path.])
      AC_MSG_NOTICE([If the package is not installed, please install it first; otherwise note])
      AC_MSG_NOTICE([that, in particular, /usr/local/lib/pkgconfig is often not in this path.])
      AC_MSG_NOTICE([Set PKG_CONFIG_PATH appropriately or use option --with-vrna=PREFIX for easy remedy.])
   fi
fi
   
if test "$VRNA_OK" = "false" -o "$VRNALIB_OK" = "false" ; then
   AC_MSG_NOTICE()
   AC_MSG_NOTICE([The Vienna RNA package (>=1.8) is required for compiling;])
   AC_MSG_NOTICE([linking (--enable-librna) requires version >= 2.1.0.])
   AC_MSG_NOTICE([It is freely available at http://www.tbi.univie.ac.at/~ivo/RNA/.])
   AC_MSG_NOTICE()
    
   exit -1;
fi


AC_OUTPUT([src/Makefile Utils/Makefile man/Makefile lib/Makefile lib/perl/Makefile lib/perl/MLocarna/Makefile Utils/mlocarna Utils/locarnate Utils/aln2fa.pl Utils/alnsel.pl Utils/pp2dot Utils/mlocarna_nnames Utils/locarna-mea Utils/reliability-profile.pl])
