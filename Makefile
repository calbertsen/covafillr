R?=R

all: doc build check install test readme

doc:
	@echo "\033[0;32mUpdating documentation\033[0;0m"
	rm -f covafillr/src/*.so
	$(R) -q -e 'devtools::document("covafillr")'

readme:
	@echo "\033[0;32mCreating README\033[0;0m"
	rm -r -f README_files
	$(R) -q -e 'rmarkdown::render("covafillr/README.Rmd",rmarkdown::md_document(variant = "gfm"))'

build: doc
	@echo "\033[0;32mBuilding package\033[0;0m"
	$(R) CMD build covafillr

check: doc build
	@echo "\033[0;32mChecking package as cran\033[0;0m"
	$(eval VERSION = $(shell Rscript -e "l<-readLines(\"covafillr/DESCRIPTION\");cat(gsub(\"Version: \",\"\",l[grepl(\"Version: \",l)]))"))
	$(R) CMD check --as-cran covafillr_${VERSION}.tar.gz

install: doc build
	@echo "\033[0;32mInstalling package\033[0;0m"
	$(eval VERSION = $(shell Rscript -e "l<-readLines(\"covafillr/DESCRIPTION\");cat(gsub(\"Version: \",\"\",l[grepl(\"Version: \",l)]))"))
	@echo "Version: ${VERSION}"
	$(R) CMD INSTALL covafillr_${VERSION}.tar.gz

test:
	@echo "\033[0;32mNo tests\033[0;0m"

clean: 
	@echo "\033[0;32mRemoving Rhistory files cleaning directory\033[0;0m"
	find . -type f -name '.Rhistory' -delete
	@echo "\033[0;32mRemoving Rcheck directory cleaning directory\033[0;0m"
	rm -f -r covafillr.Rcheck
	@echo "\033[0;32mRun package cleanup\033[0;0m"
	./covafillr/cleanup


clean_hard:
	@echo "\033[0;32mHard cleaning directory\033[0;0m"
	git clean -f -d

uninstall:
	@echo "\033[0;32mUninstalling package\033[0;0m"
	$(R) CMD REMOVE covafillr
