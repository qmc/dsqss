SUBDIRS = ../lib
AUTOMAKE_OPTIONS = subdir-objects

bin_PROGRAMS = dla_P lattgene_P


lattOBgene_P_CFLAGS =  $(INCPATH)
lattgene_P_SOURCES = lattgene.cc


dla_P_LDFLAGS = $(LDFLAGS)
dla_P_CFLAGS = \
	$(LINPATH)\
	$(INCPATH)\
	$(RNDPATH) \
	$(CFLAGS)

dla_P_LDADD = ../lib/libwrm_hrd.a

AM_CPPFLAGS = \
	-I$(top_srcdir)/src/inc \
	-I$(top_srcdir)/lib 

dla_P_SOURCES = DLA_multi.cpp \
	lattice.cpp \
	probability.cpp \
	configuration.cpp \
	quantities.cpp 
