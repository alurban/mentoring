Newtonian Tidal Distortion
==========================

In this problem, we will bring together a couple of basic things we know about orbits and make a toy model of a binary neutron star system. With nothing more than Kepler’s laws and Newton’s law of gravity, we will build an intuition for what happens as a pair of neutron stars spiral into one another, and use those insights to design an awesome orbital simulation.... of science!

Contents
--------

* `tidal_notes.{tex,pdf}`: A statement of the problem and some notes on how I solved it.
* `newtonian_tidal_forces.{py,pdf}`: A python script and its output. The script computes orbital separation and tidal forces as a function of gravitational wave frequency.
* `binary_diagram.{tex,pdf}`: A LaTeX rendering of a diagram of a neutron star binary. Pretty much does what it says on the tin.
* `matplotlibrc`: A local configuration of matplotlib settings, because I am very picky about how my plots look.
* `Makefile`: Running `make` will produce a PDF rendering of `tidal_notes.pdf`, the notes document referred to above.

Note: all the details of this problem are captured in `tidal_notes.pdf`.

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
python newtonian_tidal_forces.py
```

To produce a PDF rendering of a diagram of a binary neutron star system:

```
pdflatex binary_diagram
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
make stable_orbits
```

In the words of Captain Picard, `make` it so!
