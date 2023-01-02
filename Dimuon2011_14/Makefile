SRC_FILES=$(shell ls *.C)
FILES=$(SRC_FILES:%.C=%)

.PHONY: all $(FILES)
all: $(FILES)
$(FILES): %: %.pdf

%.pdf: %.C
	root $<+ -l -q
	@rm $*_*
