include ../Makefile.def

install:
	@if([ ! -d Python ]); then \
		git clone git\@gitlab.umich.edu:swmf_software/swmfpy Python; fi
clean:
	cd Library/src; make clean
	cd Library/test;make clean
	cd Prologs;     make clean
	rm -f include/*.mod

distclean: clean
	cd Library/src; make distclean
	cd Library/test;make distclean
	cd Prologs;     make distclean
	rm -f Library/src/mpif*.h *~ */*~

