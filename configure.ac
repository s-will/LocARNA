dnl                                               -*- Autoconf -*-
dnl Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT([LocARNA], [1.9.3beta], [will@informatik.uni-freiburg.de], [locarna],
        [http://www.bioinf.uni-freiburg.de/Software/LocARNA/])

dnl special dir for aux config files
AC_CONFIG_AUX_DIR(config)

dnl Define here to keep all version information in one place.
AC_SUBST([SO_VERSION], [1:1:0])
AC_SUBST([API_VERSION], [2.0])


dnl Package meta information
dnl
dnl vcs repository of the package
AC_SUBST([PACKAGE_VCS], ["https://github.com/s-will/LocARNA.git"])
dnl
dnl package vcs repository for browsing
AC_SUBST([PACKAGE_VCS_BROWSER], ["https://github.com/s-will/LocARNA"])
dnl
dnl main author of the package
AC_SUBST([PACKAGE_MAIN_AUTHOR], ["Sebastian Will <swill@csail.mit.edu>"])
dnl short description of the package
AC_SUBST([PACKAGE_SHORT_DESCRIPTION], ["Alignment of RNAs"])
dnl
dnl long description of the package
AC_SUBST([PACKAGE_LONG_DESCRIPTION], ["The LocARNA package consists of several \
tools for producing fast and high-quality pairwise and multiple \
alignment of RNA sequences of unknown structure. These tools build on \
the Turner free energy model of RNAs to simultaneously fold and align \
(or match) RNAs based on their sequence and structure features."])
dnl
dnl (prefix of) url of released sources (tar.gz); releases are called %{name}-%{version}.tar.gz
AC_SUBST([PACKAGE_SOURCES_URL], ["http://www.bioinf.uni-freiburg.de/Software/LocARNA/Releases"])


dnl minimally required ViennaRNA version
VRNA_VERSION=2.3.1

AC_CONFIG_MACRO_DIR([m4])

AC_CONFIG_HEADERS([src/config.h:config.h.in])

dnl include libtool
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

dnl require c++14
AX_CXX_COMPILE_STDCXX_14([noext],[mandatory])

AC_CHECK_PROGS([PERL],[perl])

dnl defines that have to be set when using the LocARNA library
AC_SUBST(LIBDEFS)

dnl Variable settings for generation of man pages via help2man --
dnl adapted from gnu hello example
dnl
dnl Are we building from vc checked-out sources, or a tarball ?
dnl This is used in "Makefile.am" to avoid re-generating the manpage
dnl when building from tarballs.
AM_CONDITIONAL([BUILD_FROM_VC], [test -d "$srcdir/.git"])
AM_CONDITIONAL([VC_CROSS_COMPILING],
               [test -d "$srcdir/.git" && test $cross_compiling = yes])
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
dnl generate doxygen documentation in Doxy
dnl

DX_HTML_FEATURE(ON)
DX_CHM_FEATURE(OFF)
DX_CHI_FEATURE(OFF)
DX_MAN_FEATURE(OFF)
DX_RTF_FEATURE(OFF)
DX_XML_FEATURE(OFF)
DX_PDF_FEATURE(OFF)
DX_PS_FEATURE(OFF)

DX_INIT_DOXYGEN($PACKAGE_NAME, ["Doc/doxygen.cfg"], ["Doxy"])
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


dnl ------------------------------
dnl Vienna RNA lib configuration
dnl
dnl check whether pkg-config is available
AC_CHECK_PROG([PKG_CONFIG_AVAILABLE],pkg-config,[true],[false])

VRNALIB_OK=false

dnl if pkg-config is available, use it to configure the Vienna RNA library

AS_IF([test "$PKG_CONFIG_AVAILABLE" = "true"],[
    dnl if with_vrna is set, then use it to set the path for pkg_config
    if test "$with_vrna" != "" ; then
        export PKG_CONFIG_PATH="$with_vrna/lib/pkgconfig:$PKG_CONFIG_PATH"
    fi
    dnl call pkg-config (macro PKG_CHECK_MODULES) to set VRNA_LIBS and VRNA_CFLAGS
    PKG_CHECK_MODULES([VRNA],[RNAlib2 >= $VRNA_VERSION],[VRNALIB_OK=true],[VRNALIB_OK=false])
])
dnl if pkg-config is not available or something else went wrong with pkg-config,
dnl then guess the configuration
dnl
AS_IF([test "$VRNALIB_OK" != "true" -a "$PKG_CONFIG_AVAILABLE" != "true" ],
[
  AC_MSG_NOTICE([Guess the Vienna RNA library configuration:])
  if test "$with_vrna" = "" ; then
    if test -e "/usr/local/include/ViennaRNA" ; then
      vrna_prefix="/usr/local"
    else
      vrna_prefix="/usr"
    fi
  else
    vrna_prefix="$with_vrna"
  fi
  VRNA_LIBS="-L$vrna_prefix/lib -lRNA"

  dnl if libgomp is available, then add it to VRNA_LIBS; otherwise,
  dnl assume (optimistically) that libRNA works without libgomp
  AC_CHECK_LIB([gomp],[main],[VRNA_LIBS="$VRNA_LIBS -lgomp"])

  VRNA_CFLAGS="-I$vrna_prefix/include"
  AC_MSG_NOTICE([Guessing of ViennaRNA flags may be avoided by installing pkg-config.])

  AC_CHECK_LIB([RNA],[vrna_fold],[VRNALIB_OK=true],[VRNALIB_OK=false])
])

# add vrna to libs
LIBS="$VRNA_LIBS $LIBS"

dnl ============================================================
dnl check for RNAlib/c++-compiler incompatibility
dnl If the ViennaRNA lib is compiled w/ C11 support, then the c++
dnl compiler has to support anonymous unions and structs.

dnl does the RNAlib require C11 support
AC_CHECK_LIB([RNA],[vrna_C11_features],[VRNA_USES_C11_FEATURES="true"],[VRNA_USES_C11_FEATURES="false"])

dnl can the C++ compiler compile anonymous unions
AC_COMPILE_IFELSE([AC_LANG_SOURCE([struct s{union{struct {float x;};struct {int i;};};};])],
                  [HAS_ANONYMOUS_STRUCT="yes"],[HAS_ANONYMOUS_STRUCT="no"])
dnl
AS_IF([test "$VRNA_USES_C11_FEATURES" = "true" && test "$HAS_ANONYMOUS_STRUCT" = "no"],
      [AC_MSG_ERROR([
      ============================================================
      Your C++ compiler does not support C11 features
      (anonymous structs) that are required by the available
      RNAlib. Please recompile the Vienna RNA library with
      --disable-c11 or use a different C++ compiler.
      ============================================================
      ])
])

dnl END check for RNAlib/c++-compiler incompatibility
dnl ============================================================

dnl disable c11 if necessary
AS_IF([test "$VRNA_USES_C11_FEATURES" = "false" ],[
    AX_APPEND_FLAG(["-DVRNA_DISABLE_C11_FEATURES"],[VRNA_CFLAGS])
])

AC_MSG_NOTICE([VRNA_LIBS=$VRNA_LIBS; VRNA_CFLAGS=$VRNA_CFLAGS])

dnl finally, set CPPFLAGS
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

dnl --------------------
dnl add linker options for POSIX threads
dnl (necessary for linking Gecode)
AX_PTHREAD
AC_MSG_NOTICE([pthread: $PTHREAD_CFLAGS, $PTHREAD_LIBS])
LIBS="$PTHREAD_LIBS $LIBS"
dnl CFLAGS="$CFLAGS $PTHREAD_CFLAGS"
LDFLAGS="$PTHREAD_CFLAGS $LDFLAGS"

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
   AX_APPEND_FLAG(["-all-static"],[LDFLAGS])
fi

dnl ----------------------------------------
dnl Static/Dynamic linking to libLocARNA
dnl
AC_MSG_CHECKING([whether to link binaries dynamically to libLocARNA])
AC_ARG_ENABLE([dynamic-liblocarna],
    AC_HELP_STRING(
        [--enable-dynamic-liblocarna],
        [link binaries dynamically to libLocARNA (def=no)]
    ),
    [dynamic_liblocarna="$enableval"],
    [dynamic_liblocarna="no"]
)
AC_MSG_RESULT([$dynamic_liblocarna])
AM_CONDITIONAL([STATIC_LIBLOCARNA],[test "$link_statically" = "no" -a "$dynamic_liblocarna" = "no"])

dnl for lgomp if all-static
AS_IF([test "$link_statically" = "yes"],
[
        AX_APPEND_FLAG(["-ldl"],[LDFLAGS])
])


dnl Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_CONST
AC_TYPE_SIZE_T


AC_CONFIG_FILES([Makefile src/Makefile src/Tests/Makefile src/Utils/Makefile Doc/Makefile Doc/man/Makefile
lib/Makefile lib/perl/Makefile lib/perl/MLocarna/Makefile
LocARNA-${API_VERSION}.pc:lib/LocARNA.pc.in
src/package_constants.hh:src/package_constants.hh.in
Packaging/debian/control Packaging/rpm/locarna.spec])


dnl ----------------------------------------
dnl Test if Vienna lib available and ok
dnl don't allow compiling otherwise

MESSAGE=""
if test "$VRNALIB_OK" = "false" ; then
   MESSAGE+="=======================================================
 "
   if test $PKG_CONFIG_AVAILABLE = "true" ; then
      if test "$with_vrna" != "" ; then
         MESSAGE+="Cannot find  pkgconfig file RNAlib2.pc for the Vienna RNA package (>= $VRNA_VERSION)
 in \'$with_vrna/lib/pkgconfig\'.
 "
      else
         MESSAGE+="Cannot find pkgconfig file RNAlib2.pc for the Vienna RNA package (>= $VRNA_VERSION)
 in the default path. Make sure that the Vienna RNA package is installed.
 If the RNA package is installed, but the pkgconfig file is still not found,
 please set PKG_CONFIG_PATH appropriately or provide --with-vrna=PREFIX.

 "
      fi
   else
         MESSAGE+="Cannot find the Vienna RNA package (>= $VRNA_VERSION) in the installation path $vrna_prefix.
 Make sure that the Vienna RNA package is installed. If the RNA package is installed,
 please specify its location --with-vrna=PREFIX.

 If possible, consider to use pkg-config as well.

 "
   fi

  MESSAGE+="LocARNA REQUIRES the Vienna RNA package (>= $VRNA_VERSION);
 It is freely available at http://www.tbi.univie.ac.at/RNA/.
"
   AC_MSG_ERROR([$MESSAGE],[-1])
fi

AC_OUTPUT
