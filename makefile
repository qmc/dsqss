
######################
ROOTDIR=$(shell pwd)
dver=1.1.18
HDIR=$(shell pwd)
pver=1.1.4

iver=8

DSQSS_INSTALL_DIR=$(ROOTDIR)/dsqss/dsqss-$(dver)
PMWA_INSTALL_DIR=$(ROOTDIR)/pmwa/pmwa-$(pver)
TOOL=$(ROOTDIR)/tool
BINDIR = $(ROOTDIR)/bin

##### DSQSS MAKER #######

EXED=$(BINDIR)/dla

##### PMWA MAKER #######

EXEP=$(BINDIR)/dla_P

##### TOOLS MAKER #######

INPGENE=inpgene-${iver}.sh

########################

LATS=\
	$(BINDIR)/lattgene_P \
	$(BINDIR)/lattgene \

ALGS=\
	$(BINDIR)/hamgen_H \

HAMS=\
	$(BINDIR)/dla_alg \

########################

all:
	$(MAKE) -C $(DSQSS_INSTALL_DIR)
	$(MAKE) -C $(PMWA_INSTALL_DIR)
	$(MAKE) -C $(DSQSS_INSTALL_DIR) install
	$(MAKE) -C $(PMWA_INSTALL_DIR) install
	cp $(TOOL)/$(INPGENE) $(BINDIR)/inpgene
	cp $(TOOL)/run_merge.sh $(BINDIR)/runmerge
	cp $(TOOL)/extrap.pl $(BINDIR)/extrap.pl
	cp -r  $(TOOL)/plot_20150806 $(BINDIR)/plot

install:
	$(MAKE) -C $(DSQSS_INSTALL_DIR) install
	$(MAKE) -C $(PMWA_INSTALL_DIR) install
clean:
	$(MAKE) -C $(DSQSS_INSTALL_DIR) clean
	$(MAKE) -C $(PMWA_INSTALL_DIR) clean
	$(RM) $(EXED) $(EXEP) $(LATS) $(ALGS) $(HAMS)

check:
	$(MAKE) -C $(DSQSS_INSTALL_DIR) check

lattgene:
	$(MAKE) -C $(DSQSS_INSTALL_DIR) lattgene

dlg_alg:
	$(MAKE) -C $(DSQSS_INSTALL_DIR) dla_alg
dist:
	$(MAKE) -C $(DSQSS_INSTALL_DIR) dist
distx:
	$(MAKE) -C $(DSQSS_INSTALL_DIR) distx
distx2:
	$(MAKE) -C $(DSQSS_INSTALL_DIR) distx2
sphinx:
	$(MAKE) -C $(DSQSS_INSTALL_DIR) sphinx
doxygen:
	$(MAKE) -C $(DSQSS_INSTALL_DIR) doxygen
