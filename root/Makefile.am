#
# $Id: Makefile.am,v 1.14 2008/02/12 16:00:06 barbui Exp $
#
AUTOMAKE_OPTIONS        = foreign

#CLEANFILES              = *~ core G__autoLinkDef.h *Cint.* *.root \
#                          tmp.C .rootrc

MAINTAINERCLEANFILES    = Makefile.in Makefile configure config.* \
                          missing mkinstalldirs libtool ltmain.sh \
                          ltconfig.sh install-sh depcomp \
			  AUTHORS COPYING INSTALL NEWS README ChangeLog\
			  aclocal.m4 

%Cint.cxx %Cint.h:Include.h LinkDef.h
	$(ROOTCINT) -f $*Cint.cxx -c $(AM_CPPFLAGS) $(INCLUDES) $(DEFS) $^

DICTSOURCES=GeminiCint.cxx GeminiCint.h LinkDef.h

lib_LTLIBRARIES		= libGemini.la

noinst_HEADERS      	= Include.h LinkDef.h

AM_CPPFLAGS		= -I$(ROOTINCDIR) -I$(HOME)/include

libGemini_la_SOURCES =  Nucleus.cpp			\
			Mass.cpp			\
			Chart.cpp 			\
			Yrast.cpp			\
			TlArray.cpp			\
			LevelDensity.cpp		\
			Angle.cpp			\
			Nuclide.cpp			\
			LightP.cpp			\
			Evap.cpp			\
			AngleDist.cpp			\
			Random.cpp			\
                        Scission.cpp                    \
                        TlBarDist.cpp                   \
                        Fus.cpp                         \
                        Weight.cpp                      \
                        SigCharged.cpp                  \
                        SigBarDist.cpp                  \
                        Gdr.cpp                         \
			$(DICTSOURCES)

pkginclude_HEADERS =	CNucleus.h

# Change Here Too!!
libGemini_la_LDFLAGS	= -L$(HOME)/lib -module -avoid-version \
			  -R @ROOTLIBDIR@

bin_PROGRAMS		= testDecay

testDecay_LDFLAGS 	= $(ROOTLIBS) -R $(ROOTRPATH) -L$(ROOTLIBDIR) libGemini.la

testDecay_SOURCES 	= testDecay.cpp


