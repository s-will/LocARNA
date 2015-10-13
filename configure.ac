#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT([LocARNA], [1.8.5b], [will@informatik.uni-freiburg.de], [locarna],
        [http://www.bioinf.uni-freiburg.de/Software/LocARNA/])

# Define here to keep all version information in one place.
AC_SUBST([SO_VERSION], [1:0:0])
AC_SUBST([API_VERSION], [1.3])

AC_CONFIG_MACRO_DIR([m4])

AC_CONFIG_HEADERS([src/config.h:config.h.in])

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

AC_CHECK_PROGS([PERL],[perl])
AC_CHECK_PROG([HELP2MAN],[help2man],[help2man],[false])
AC_CHECK_PROG([POD2MAN],[pod2man],[pod2man],[false])

## ----------------------------------------
# doxygen documentation
#

DX_HTML_FEATURE(ON)
DX_CHM_FEATURE(OFF)
DX_CHI_FEATURE(OFF)
DX_MAN_FEATURE(OFF)
DX_RTF_FEATURE(OFF)
DX_XML_FEATURE(OFF)
DX_PDF_FEATURE(OFF)
DX_PS_FEATURE(OFF)

DX_INIT_DOXYGEN($PACKAGE_NAME, ["doxygen.cfg"], ["Doc"])
## ----------------------------------------



## ----------------------------------------
# for Vienna RNA
AC_CHECK_FUNCS([strdup]) #defines HAVE_STRDUP if strdup is available

# Vienna RNA package library path support, if not installed in standard directory/path
VRNA_OK=true

AC_ARG_WITH([vrna],
    AC_HELP_STRING(
        [--with-vrna=PREFIX],
        [Installation directory of the Vienna RNA library]
      )
)

ENABLE_LIBRNA=yes ## require rna lib linking

AC_SUBST(LDFLAGS,[""])
AC_SUBST(CPPFLAGS,[""])
AC_SUBST(LIBS,[""])
AC_SUBST(LIBDEFS,[""]) # defines that have to be set when using the LocARNA library


## ------------------------------
## Vienna RNA lib configuration
##
## check whether pkg-config is available 
AC_CHECK_PROG([PKG_CONFIG_AVAILABLE],pkg-config,[true],[false])

VRNALIB_OK=false

##
## if pkg-config is available, use it to configure the Vienna RNA library
AS_IF([test "$PKG_CONFIG_AVAILABLE" = "true"],[
    # set PKG_CONFIG_PATH
    if test "$with_vrna" != "" ; then
        export PKG_CONFIG_PATH="$with_vrna/lib/pkgconfig:$PKG_CONFIG_PATH"
        if test ! -e "$with_vrna/lib/pkgconfig/RNAlib2.pc" ; then
            AC_MSG_NOTICE([Pkg-config of Vienna RNA cannot be found.])
            AC_MSG_NOTICE([Fall back to guessing.])
            VRNA_GUESS_FALLBACK="true"
        fi
    fi
    #
    AS_IF([test "$VRNA_GUESS_FALLBACK" != "true"],[
        PKG_CHECK_MODULES([VRNA],[RNAlib2 >= 2.1.1],[VRNALIB_OK=true],[VRNALIB_OK=false])
        if test "$with_vrna" != "" ; then
            VRNA_CFLAGS="$VRNA_CFLAGS -I$with_vrna/include"
        fi
    ])
])
## if pkg-config is not available or something else went wrong with pkg-config,
## then guess the configuration
##
AS_IF([test "$VRNALIB_OK" != "true"], 
[
  AC_MSG_NOTICE([Guess the Vienna RNA library configuration:])
  if test "$with_vrna" = "" ; then
    vrna_prefix="/usr/local"
  else
    vrna_prefix="$with_vrna"
  fi
  VRNA_LIBS="-L$vrna_prefix/lib -lRNA"

  ## if libgomp is available, then add it to VRNA_LIBS; otherwise,
  ## assume (optimistically) that libRNA works without libgomp 
  AC_CHECK_LIB([gomp],[main],[VRNA_LIBS="$VRNA_LIBS -lgomp"])
  
  VRNA_CFLAGS="-I$vrna_prefix/include"
  AC_MSG_NOTICE([VRNA_LIBS=$VRNA_LIBS; VRNA_CFLAGS=$VRNA_CFLAGS])
  AC_MSG_NOTICE([Guessing may be avoided by installing pkg-config.])
  VRNALIB_OK=true
])
##
## finally, set LDFLAGS and CPPFLAGS
LDFLAGS="$VRNA_LIBS $LDFLAGS"
CPPFLAGS="$VRNA_CFLAGS $CPPFLAGS"

# ------------------------------
# handle open mp support
#  this is necessary for linking against the vienna lib 2.1.x
#AC_OPENMP
#AC_SUBST(OPENMP_CXXFLAGS)
#
#CPPFLAGS="$OPENMP_CXXFLAGS $CPPFLAGS"
#LDFLAGS="$OPENMP_CXXFLAGS $LDFLAGS"

# Checks for header files.
AC_HEADER_STDC

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
if test "$debug_it" = "yes"; then
    AC_DEFINE([DEBUG],1,[Turn on debugging])
    CPPFLAGS="$CPPFLAGS -g -Wall"
else
    AC_DEFINE([NDEBUG],1,[Turn off debugging])
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
if test "$link_statically" = "yes"; then
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
if test "$use_large_pf" = "yes"; then
    LIBDEFS="$LIBDEFS -DVERY_LARGE_PF"
    AC_DEFINE([VERY_LARGE_PF],1,[Use long double type for partition functions])
fi

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_CONST
AC_TYPE_SIZE_T


AC_CONFIG_FILES([Makefile src/Makefile Utils/Makefile man/Makefile lib/Makefile lib/perl/Makefile lib/perl/MLocarna/Makefile Utils/exploc_p LocARNA-${API_VERSION}.pc:LocARNA.pc.in])


### ----------------------------------------
### Test if Vienna lib available and ok
## don't allow compiling otherwise 
#

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
   
if test "$VRNALIB_OK" = "false" ; then
   AC_MSG_NOTICE()
   AC_MSG_NOTICE([The Vienna RNA package (>=2.1.1) is required;])
   AC_MSG_NOTICE([It is freely available at http://www.tbi.univie.ac.at/RNA/.])
   AC_MSG_NOTICE()
    
   exit -1;
fi

AC_OUTPUT
