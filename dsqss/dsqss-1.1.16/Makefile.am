SUBDIRS = src 
EXTRA_DIST = \
runConfigure \
./tests/* \
./tool/* \
./samp/* \
./include/* \
./doc/* 

#dist_data = `date + %Y%m%d%H%M`
ams_date = `date +%Y%m%d`
ams_date2 = `date +%Y%m%d%H%M`
distx:
	$(MAKE) dist
	mv $(distdir).tar.bz2 $(distdir)-$(ams_date).tar.bz2

distx2:
	$(MAKE) dist
	mv $(distdir).tar.bz2 $(distdir)-$(ams_date2).tar.bz2

sub: lattgene dla_alg hamgen_H

lattgene:
	cd src; $(CXX) $(CPPFLAGS) -pedantic-errors -o lattgene lattgene.cc

dla_alg: dla_alg.o
	cd src; $(CXX) $(CPPFLAGS) -pedantic-errors -o $@ dla_alg.o 

dla_alg.o:
	cd src; $(CXX) $(CPPFLAGS) -pedantic-errors -c dla_alg.cc io.h xml.h array.h

#hamgen_H: exact_H
#	cd src; ln -sf exact_H hamgen_H

hamgen_H:
	cd src; $(CXX) $(CPPFLAGS) $(CXXFLAGS) -pedantic-errors -o exact_H exact_H.cc $(LIBS)

doxygen:
	cd src; doxygen doxyfile.in

sphinx:
	cd doc/users_manual ;\
	make html;
	cd doc/users_manual/build/html ;\
	cat index.html | sed -n '1,4p' >> aaa ;\
	cat ../../source/counter.html >> aaa ;\
	more +6 index.html >> aaa ;\
	mv aaa index.html

TESTS = tests/tests.sh
