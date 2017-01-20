BYPROD = *.aux *.log *.out *.dvi *.blg *.bbl *.xml *.bcf *.synctex.gz *Notes.bib

all: tidal_notes stable_orbit

%: %/*.pdf %.tex
	pdflatex -draftmode $(@:.pdf=) && \
	pdflatex -draftmode $(@:.pdf=) && \
	pdflatex -draftmode $(@:.pdf=) && \
	pdflatex $(@:.pdf=)


clean:
	rm -f $(BYPROD)
