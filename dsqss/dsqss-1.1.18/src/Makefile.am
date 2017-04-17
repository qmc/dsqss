bin_PROGRAMS = dla dla_alg cfgene sfgene lattgene_C lattgene_T hamgen_H hamgen_B
#bin_PROGRAMS =  exact_H
dla_SOURCES = \
  random.cc \
  dla.cc 

dla_alg_SOURCES = dla_alg.cc
cfgene_SOURCES = cfgene.cc
sfgene_SOURCES = sfgene.cc
lattgene_C_SOURCES = lattgene_C.cc
lattgene_T_SOURCES = lattgene_T.cc
hamgen_H_SOURCES = exact_H.cc
hamgen_B_SOURCES = exact_B.cc

#exact_H_LDFLAGS= -I/home/kota/local/lib/cpplapack/include
#exact_H_CXXFLAGS= -I/home/kota/local/lib/cpplapack/include
#exact_H_INCLUDES= -I/home/kota/local/lib/cpplapack/include

#install-exec-local:
#	ln -sf exact_H ${bindir}/hamgen_H

EXTRA_DIST=\
graphic.cc \
algorithm.hpp \
dla.hpp \
link.hpp \
objects.hpp \
chainjob.hpp \
lattice.hpp \
measure.hpp \
parameter.hpp \
dsqss_read.hpp \
array.h \
io.h \
measure_specific.h \
name.h \
random.h \
xml.h \
dla_alg.h \
canonical.h \
matrix.h \
spin_H.h \
boson_B.h \
measure_specific.cc 
