Tidal Distortion
================

In this problem, we will bring together a couple of basic things we know about orbits and make a toy model of a binary neutron star system. With nothing more than Kepler’s laws and Newton’s law of gravity, we will build an intuition for what happens as a pair of neutron stars spiral into one another, and use those insights to design an awesome orbital simulation.... of science!

Contents
--------

* `tidal_notes.{tex,pdf}`: The first problem set on Newtonian tidal forces, and some notes on how I solved it.
* `keplerian_orbit.{tex,pdf}`: The second problem set on stable Keplerian orbits, and some notes on how I solved it.
* `scripts/`: A directory containing Python scripts related to these notes, e.g. a script that computes orbital separation and tidal forces as a function of gravitational wave frequency.
* `scripts/matplotlibrc`: A local configuration of matplotlib settings, because I am very picky about how my plots look.
* `Makefile`: Running `make` will produce a PDF rendering of the notes documents referred to above.

Note: all the details of this problem are captured in `*.pdf` documents in this directory. In sequntial order, they go:

1. `tidal_notes.pdf`
2. `keplerian_orbit.pdf`

Software Requirements
---------------------

* [LaTeX 2e](https://www.latex-project.org/get/)
* pdflatex
* [Python 2.x](https://www.python.org)
* [NumPy](http://www.numpy.org)
* [Matplotlib](http://matplotlib.org)

Instructions
------------

To produce a plot of the orbital separation and tidal forces as a function of gravitational wave frequency:

```
python scripts/newtonian_tidal_forces.py
```

To produce a PDF rendering of a diagram of a binary neutron star system:

```
pdflatex scripts/binary_diagram
```

To produce a PDF rendering of my complete notes on this problem, with suggested solution and thinking points:

```
pdflatex -draftmode tidal_notes && \
pdflatex -draftmode tidal_notes && \
pdflatex tidal_notes
```

and the same for notes on stable orbits. Alternatively, you can use the Makefile, which `make`s life much simpler:

```
make tidal_notes
make keplerian_orbits
```

To re-generate drafts of all notes:

```
make
```

In the words of Captain Picard, `make` it so!
