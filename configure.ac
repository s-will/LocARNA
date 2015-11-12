dnl                                               -*- Autoconf -*-
dnl Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT([LocARNA], [1.8.7], [will@informatik.uni-freiburg.de], [locarna],
        [http://www.bioinf.uni-freiburg.de/Software/LocARNA/])

dnl Define here to keep all version information in one place.
AC_SUBST([SO_VERSION], [1:0:0])
AC_SUBST([API_VERSION], [1.4])

AC_CONFIG_MACRO_DIR([m4])

AC_CONFIG_HEADERS([src/config.h:config.h.in])

dnl for libtool
LT_INIT
dnl deprecated; left for backward compatibility
AC_PROG_LIBTOOL

AM_INIT_AUTOMAKE

dnl Checks for programs.
AC_PROG_CXX
AC_PROG_CXXCPP
AC_PROG_CC
AC_PROG_CPP

AC_LANG_CPLUSPLUS

AC_CHECK_PROGS([PERL],[perl])

dnl Variable settings for generation of man pages via help2man --
dnl adapted from gnu hello example
dnl 
dnl Are we building from vc checked-out sources, or a tarball ?
dnl This is used in "Makefile.am" to avoid re-generating the manpage
dnl when building from tarballs.
AM_CONDITIONAL([BUILD_FROM_VC], [test -d "$srcdir/.hg"])
AM_CONDITIONAL([VC_CROSS_COMPILING],
               [test -d "$srcdir/.hg" && test $cross_compiling = yes])
dnl
dnl GNU help2man creates man pages from --help output; in many cases, this
dnl is sufficient, and obviates the need to maintain man pages separately.
dnl However, this means invoking executables, which we generally cannot do
dnl when cross-compiling, so we test to avoid that (the variable
dnl "cross_compiling" is set by AC_PROG_CC).
AS_IF([test $cross_compiling = no], [
  AM_MISSING_PROG([HELP2MAN], [help2man])
], [
  HELP2MAN=:
])

AC_CHECK_PROG([POD2MAN],[pod2man],[pod2man],[false])

dnl ----------------------------------------
dnl doxygen documentation
dnl

DX_HTML_FEATURE(ON)
DX_CHM_FEATURE(OFF)
DX_CHI_FEATURE(OFF)
DX_MAN_FEATURE(OFF)
DX_RTF_FEATURE(OFF)
DX_XML_FEATURE(OFF)
DX_PDF_FEATURE(OFF)
DX_PS_FEATURE(OFF)

DX_INIT_DOXYGEN($PACKAGE_NAME, ["doxygen.cfg"], ["Doc"])
dnl ----------------------------------------



dnl ----------------------------------------
dnl for Vienna RNA
dnl     defines HAVE_STRDUP if strdup is available
AC_CHECK_FUNCS([strdup]) 

dnl Vienna RNA package library path support, if not installed in standard directory/path
VRNA_OK=true

AC_ARG_WITH([vrna],
    AC_HELP_STRING(
        [--with-vrna=PREFIX],
        [Installation directory of the Vienna RNA library]
      )
)

dnl require rna lib linking
ENABLE_LIBRNA=yes

AC_SUBST(LDFLAGS,[""])
AC_SUBST(CPPFLAGS,[""])
AC_SUBST(LIBS,[""])
dnl defines that have to be set when using the LocARNA library
AC_SUBST(LIBDEFS,[""]) 


dnl ------------------------------
dnl Vienna RNA lib configuration
dnl
dnl check whether pkg-config is available 
AC_CHECK_PROG([PKG_CONFIG_AVAILABLE],pkg-config,[true],[false])

VRNALIB_OK=false

dnl if pkg-config is available, use it to configure the Vienna RNA library

AS_IF([test "$PKG_CONFIG_AVAILABLE" = "true"],[
    if test "$with_vrna" != "" ; then
        export PKG_CONFIG_PATH="$with_vrna/lib/pkgconfig:$PKG_CONFIG_PATH"
        if test ! -e "$with_vrna/lib/pkgconfig/RNAlib2.pc" ; then
            AC_MSG_NOTICE([Pkg-config of Vienna RNA cannot be found.])
            AC_MSG_NOTICE([Fall back to guessing.])
            VRNA_GUESS_FALLBACK="true"
        fi
    fi
    AS_IF([test "$VRNA_GUESS_FALLBACK" != "true"],[
        PKG_CHECK_MODULES([VRNA],[RNAlib2 >= 2.1.1],[VRNALIB_OK=true],[VRNALIB_OK=false])
        if test "$with_vrna" != "" ; then
            VRNA_CFLAGS="$VRNA_CFLAGS -I$with_vrna/include"
        fi
    ])
])
dnl if pkg-config is not available or something else went wrong with pkg-config,
dnl then guess the configuration
dnl 
AS_IF([test "$VRNALIB_OK" != "true"], 
[
  AC_MSG_NOTICE([Guess the Vienna RNA library configuration:])
  if test "$with_vrna" = "" ; then
    vrna_prefix="/usr/local"
  else
    vrna_prefix="$with_vrna"
  fi
  VRNA_LIBS="-L$vrna_prefix/lib -lRNA"

  dnl if libgomp is available, then add it to VRNA_LIBS; otherwise,
  dnl assume (optimistically) that libRNA works without libgomp 
  AC_CHECK_LIB([gomp],[main],[VRNA_LIBS="$VRNA_LIBS -lgomp"])
  
  VRNA_CFLAGS="-I$vrna_prefix/include"
  AC_MSG_NOTICE([VRNA_LIBS=$VRNA_LIBS; VRNA_CFLAGS=$VRNA_CFLAGS])
  AC_MSG_NOTICE([Guessing may be avoided by installing pkg-config.])
  VRNALIB_OK=true
])

dnl finally, set LDFLAGS and CPPFLAGS
LDFLAGS="$VRNA_LIBS $LDFLAGS"
CPPFLAGS="$VRNA_CFLAGS $CPPFLAGS"

dnl Checks for header files.
AC_HEADER_STDC

dnl --------------------
dnl Debug option
dnl
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

dnl  CPPFLAGS="$CPPFLAGS -Wno-deprecated"



dnl ----------------------------------------
dnl Static linking
dnl
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


dnl flag for use of long double partition function type
dnl attention: name clash with Vienna package macro LARGE_PF ==> define VERY_LARGE_PF
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

dnl Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_CONST
AC_TYPE_SIZE_T


AC_CONFIG_FILES([Makefile src/Makefile Utils/Makefile man/Makefile lib/Makefile lib/perl/Makefile lib/perl/MLocarna/Makefile LocARNA-${API_VERSION}.pc:LocARNA.pc.in])


dnl ----------------------------------------
dnl Test if Vienna lib available and ok
dnl don't allow compiling otherwise 

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
