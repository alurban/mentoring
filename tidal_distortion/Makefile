FIGS = binary_diagram.pdf newtonian_tidal_forces.pdf
BYPROD = *.aux *.log *.out *.dvi *.blg *.bbl *.xml *.bcf *.synctex.gz *Notes.bib

all: tidal_notes.pdf

%: $(FIGS) %.tex
	pdflatex -draftmode $(@:.pdf=) && \
	pdflatex -draftmode $(@:.pdf=) && \
	pdflatex -draftmode $(@:.pdf=) && \
	pdflatex $(@:.pdf=)


clean:
	rm -f $(BYPROD)
