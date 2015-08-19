dnl 
dnl $Id: acinclude.m4,v 1.1.1.1 2007/08/28 15:53:59 psahu Exp $
dnl
dnl Macro to turn on/off debugging

AC_DEFUN([GEMINI_DEBUG], [
  dnl AC_REQUIRE(AC_PROG_CC)
  dnl AC_REQUIRE(AC_PROG_CXX)
	
  AC_MSG_CHECKING(wether to make debug objects) 

  AC_ARG_ENABLE(debug, 
    [  --enable-debug          Enable debugging symbols in objects],
      user_debug=$enableval, user_debug=yes )

  if test "x$user_debug" = "xno" ; then 
    # disable debugging symbols 
    CFLAGS=`echo $CFLAGS | sed 's,-g,,'`
    CXXFLAGS=`echo $CXXFLAGS | sed 's,-g,,'`
  fi

  AC_MSG_RESULT($user_debug 'CFLAGS=$CFLAGS')
])

AC_DEFUN([GEMINI_OPTIMIZE], [
  dnl AC_REQUIRE(AC_PROG_CC)
  dnl AC_REQUIRE(AC_PROG_CXX)
	
  AC_ARG_ENABLE(optimization, 
    [  --enable-optimization   Enable optimization of objects],
      user_optim=$enableval, user_optim=yes )

  AC_MSG_CHECKING(for optimiztion level) 
    
  if test "x$user_optim" = "xno" ; then 
    # disable debugging symbols 
    CFLAGS=`echo $CFLAGS | sed 's,-O2,,'`
    CXXFLAGS=`echo $CXXFLAGS | sed 's,-O2,,'`
  elif test "x$user_optim" = "xyes" ; then 
    # enable optimisation, do nothing
    case $host in 
    *-*-linux-*) 
      if uname -a | grep -q SMP ; then
	CFLAGS=`echo $CFLAGS | sed 's,-O2,-O1,'`
	CXXFLAGS=`echo $CXXFLAGS | sed 's,-O2,-O1,'`
      fi
      ;;
    esac 
  else 
    # Custom optimisation level 
    CFLAGS=`echo $CFLAGS | sed "s,-O2,-O$user_optim,"`
    CXXFLAGS=`echo $CXXFLAGS | sed "s,-O2,-O$user_optim,"`
  fi
 
  AC_MSG_RESULT($user_optim 'CFLAGS=$CFLAGS')
])

		
AC_DEFUN([ROOT_PATH],
[
    AC_ARG_WITH(rootsys,
    [  --with-rootsys          top of the ROOT installation directory],
	user_rootsys=$withval,
	user_rootsys="none")
    if test ! x"$user_rootsys" = xnone; then
	rootbin="$user_rootsys/bin"
    elif test ! x"$ROOTSYS" = x ; then 
	rootbin="$ROOTSYS/bin"
    else 
	rootbin=$PATH
    fi
    AC_PATH_PROG(ROOTCONF, root-config , no, $rootbin)
    if test x"$ROOTCONF" = "xno" ; then
	AC_MSG_ERROR([ROOT config script not found!])
    fi

    AC_PATH_PROG(ROOTEXEC, root, no, $rootbin) 
	
    changequote(<<, >>)dnl
    ROOTCINT=`dirname $ROOTCONF`
    ROOTCINT="${ROOTCINT}/rootcint"
    ROOTLIBDIR=`$ROOTCONF --libdir`
    ROOTINCDIR=`$ROOTCONF --incdir`
dnl ROOTLIBDIR=`$ROOTCONF --prefix`/lib
dnl ROOTINCDIR=`$ROOTCONF --prefix`/include/root
    ROOTCFLAGS=`$ROOTCONF --noauxcflags --cflags` 
dnl ROOTCFLAGS=`$ROOTCONF --cflags`
    ROOTLDFLAGS=`$ROOTCONF --libs | sed -n -e's,-l[^[:space:]]*,,gp'` 
dnl ROOTLIBS=`$ROOTCONF --noauxlibs --noldflags --libs`
    ROOTLIBS=`$ROOTCONF --libs`
dnl ROOTGLIBS=`$ROOTCONF --noauxlibs --noldflags --glibs`
    ROOTGLIBS=`$ROOTCONF --glibs`
    ROOTAUXCFLAGS=`$ROOTCONF --auxcflags`
dnl ROOTAUXCFLAGS=`$ROOTCONF --cflags`
    ROOTAUXLIBS=`$ROOTCONF --auxlibs`
dnl ROOTAUXLIBS=`$ROOTCONF --libs`
dnl ROOTMACRODIR=`$ROOTCONF --prefix`/share/root/macros
    ROOTRPATH=$ROOTLIBDIR
    changequote([, ])dnl

    AC_SUBST(ROOTCINT)
    AC_SUBST(ROOTLIBDIR)
    AC_SUBST(ROOTINCDIR)
    AC_SUBST(ROOTCFLAGS)
    AC_SUBST(ROOTLDFLAGS)
    AC_SUBST(ROOTLIBS)
    AC_SUBST(ROOTGLIBS) 
    AC_SUBST(ROOTAUXLIBS)
    AC_SUBST(ROOTAUXCFLAGS)
    AC_SUBST(ROOTRPATH)
    AC_SUBST(ROOTMACRODIR)
])
