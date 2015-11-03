
######################
ROOTDIR=$(shell pwd)
dver=1.1.17
HDIR=$(shell pwd)
pver=1.1.2

iver=7

DSQSS_INSTALL_DIR=$(ROOTDIR)/dsqss/dsqss-$(dver)
PMWA_INSTALL_DIR=$(ROOTDIR)/pmwa/pmwa-$(pver)
TOOL=$(ROOTDIR)/tool
BINDIR = $(ROOTDIR)/bin

##### DSQSS MAKER #######

EXED=$(BINDIR)/dla

##### PMWA MAKER #######

EXEP=$(BINDIR)/pmwa

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
	cp -r  $(TOOL)/plot_20150806 $(BINDIR)/plot
clean:
	$(MAKE) -C $(DSQSS_INSTALL_DIR) clean
	$(MAKE) -C $(PMWA_INSTALL_DIR) clean
	$(RM) $(EXED) $(EXEP) $(LATS) $(ALGS) $(HAMS)
