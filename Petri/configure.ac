AC_INIT([PetriSpot],[0.1],[Yann.Thierry-mieg@lip6.fr])
# We use automake to build makefiles

AM_INIT_AUTOMAKE([-Wall -Werror foreign subdir-objects])


# option to indicate libddd folder
AC_ARG_WITH([libddd],
            [AS_HELP_STRING([--with-libddd=/root/of/libddd],[specify the path of libddd])])


AM_CONDITIONAL([WITH_LIBDDD_PATH], [test "x${with_libddd}" != x])
if test "x${with_libddd}" != x; then
    AC_SUBST([LIBDDD_ROOT],[${with_libddd}])
fi

# link time opt, set AR and RANLIB before LT_INIT call, unless configure overrides them.
# also Options for disabling LTO on Darwin since binutils sucks big time on that platform.
nolto=false

AC_ARG_ENABLE(	[nolto],
        [AS_HELP_STRING([--enable-nolto],[avoid using LTO flags, mostly for Darwin])],
        [  case "${enable_nolto}" in
            yes) nolto=true
                ;;
            no)
                ;;
            *)  AC_MSG_ERROR(Bad value ${enableval})
                ;;
            esac])

if test "x$nolto" = "xfalse"; then
   CXXFLAGS="-flto $CXXFLAGS"
   test -z "$AR" && AR=gcc-ar
   test -z "$RANLIB" && RANLIB=gcc-ranlib
fi


# option to indicate libspot folder
AC_ARG_WITH([libspot],
            [AS_HELP_STRING([--with-libspot=/root/of/libspot],[specify the path of Spot])])


AM_CONDITIONAL([WITH_LIBSPOT_PATH], [test "x${with_libspot}" != x])
if test "x${with_libspot}" != x; then
    AC_SUBST([LIBSPOT_ROOT],[${with_libspot}])
fi


# option to indicate libexpat folder
AC_ARG_WITH([libexpat],
            [AS_HELP_STRING([--with-libexpat=/root/of/libexpat],[specify the path of libexpat])])


AM_CONDITIONAL([WITH_LIBEXPAT_PATH], [test "x${with_libexpat}" != x])
if test "x${with_libexpat}" != x; then
    AC_SUBST([LIBEXPAT_ROOT],[${with_libexpat}])
fi

# option to indicate the path to the antlr C runtime
AC_ARG_WITH([antlrc],
            [AS_HELP_STRING([--with-antlrc=/root/of/antlrc],[specify the path of antlr C runtime])])

AM_CONDITIONAL([WITH_ANTLRC_PATH], [test "x${with_antlrc}" != x])
if test "x${with_antlrc}" != x; then
    AC_SUBST([ANTLRC_ROOT],[${with_antlrc}])
fi



# AC_CHECK_LIB(dl, dlopen, ,AC_MSG_WARN(No operational dynamic link libdl found. This is a real problem.))

#Options for disabling google hash
stdhash=false

AC_ARG_ENABLE(	[stdhash],
        [AS_HELP_STRING([--enable-stdhash],[revert std lib hash map(default uses google sparse hash)])],
        [  case "${enable_stdhash}" in
            yes) stdhash=true
                ;;
            no)
                ;;
            *)  AC_MSG_ERROR(Bad value ${enableval})
                ;;
            esac])

if test "x$stdhash" = "xtrue"; then
    AC_DEFINE([USE_STD_HASH],1,[Define if you want to use std hash instead of google sparse hash.])
fi

CPPFLAGS="-Wno-unused-local-typedefs $CPPFLAGS"
AC_SUBST([STATICFLAGS], ["-all-static -static-libgcc -static-libstdc++"])

AC_ARG_ENABLE([mingw-native],
				[AS_HELP_STRING([--enable-mingw-native],[build native win32 binaries for redistribution])],
				[ case "${enable_mingw_native}" in
					 yes) CXXFLAGS="-DPSAPI_VERSION=1  $CXXFLAGS"
					 LDFLAGS=" -Wl,-static,--stack,16777216 $LDFLAGS"
					 LIBS="-lPSAPI -lstdc++ -lpthread $LIBS"
					 ;;
					 no)
					 ;;
					 *) AC_MSG_ERROR(Bad value ${enableval})
					 ;;
				esac])



# programs used to build
# we make static libs (use ranlib)
AC_PROG_RANLIB

# we make shared libs (use libtool)
LT_INIT([disable-shared])

# no dynamic libs unless --enable-shared is passed to configure
AC_DISABLE_SHARED
# build static libs by default
AC_ENABLE_STATIC


_LT_SET_OPTION([LT_INIT],[dlopen])
m4_warn([obsolete],[AC_LIBTOOL_DLOPEN: Remove this warning and the call to _LT_SET_OPTION when you
put the 'dlopen' option into LT_INIT's first parameter.])

# we compile C++
# set Cxx compilation flags to either the value passed in
# ./configure --cxxflags=
# or none to change behavior from default cxx flags generated by autoconf
# default : -02 -g
test -z "$CFLAGS" && CFLAGS=
test -z "$CXXFLAGS" && CXXFLAGS=

# we compile C++20
CXXFLAGS="-std=c++20 $CXXFLAGS"

#compile in C++
AC_PROG_CXX

AC_CHECK_TOOL([STRIP],[strip])

# automake :
# build these Makefile using the Makefile.in files in each directory
AC_CONFIG_FILES([
  Makefile
  src/Makefile
])

# Do it !! (mandatory)
AC_OUTPUT
