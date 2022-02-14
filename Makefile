RMD := $(wildcard *.Rmd)
HTML := $(RMD:.Rmd=.html)
MD := $(RMD:.Rmd=.md)

all: $(HTML) $(MD)

clean:
	rm -f $(MD) $(HTML)

%.html: %.Rmd %.md
	Rscript -e "rmarkdown::render('$<', output_format='html_document', output_file='$@');"

%.md: %.Rmd
	Rscript -e "rmarkdown::render('$<', output_format='md_document', output_file='$@');"
