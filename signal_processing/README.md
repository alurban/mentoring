Gravitational Wave Astronomy with LIGO and Virgo
================================================

In this directory you will find notes (and eventually, a suite of tools we develop) for introductory material in gravitational wave astronomy.
Our methods will probably focus exclusively on the ground-based detectors LIGO and Virgo. We'll begin by looking at some back-of-the-envelope
astrophysics, which will hopefully motivate our exploration of how the detectors work and how we measure properties of astrophysical sources
through their observed gravitational wave signals.

Contents
--------

* `intro.{tex,pdf}`: The first problem set introducing us to gravitational wave astronomy.
* `Makefile`: Running `make` will produce a PDF rendering of the notes documents referred to above.

Note: all the details of this problem are captured in `*.pdf` documents in this directory. In sequntial order, they go:

1. `intro.pdf`

Software Requirements
---------------------

* [LaTeX 2e](https://www.latex-project.org/get/)
* pdflatex
* [Python 2.x](https://www.python.org)
* [NumPy](http://www.numpy.org)
* [Matplotlib](http://matplotlib.org)
* [Image Magick](http://www.imagemagick.org/script/index.php), a suite of open source tools I used to turn time-sampled simulation data into gifs illustrating orbits
* [FFmpeg](https://ffmpeg.org), a video format conversion architecture that can help produce much more memory-savvy video files

Instructions
------------

To produce a PDF rendering of my notes, with suggested solution and thinking points:

```
pdflatex -draftmode intro && \
pdflatex -draftmode intro && \
pdflatex intro
```

Alternatively, you can use the Makefile, which `make`s life much simpler:

```
make intro
```

To re-generate drafts of all notes:

```
make
```
