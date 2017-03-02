Tidal Distortion
================

In this problem, we will bring together a couple of basic things we know about orbits and make a toy model of a binary neutron star system. With nothing more than Kepler’s laws and Newton’s law of gravity, we will build an intuition for what happens as a pair of neutron stars spiral into one another, and use those insights to design an awesome orbital simulation.... of science!

Contents
--------

* `intro.{tex,pdf}`: The first problem set on Newtonian tidal forces, and some notes on how I solved it.
* `keplerian_orbit.{tex,pdf}`: The second problem set on stable Keplerian orbits, and some notes on how I solved it.
* `relativistic_orbit.{tex,pdf}`: The third problem set that introduces a more efficient numerical integration scheme called RK4 (RK for Runge-Kutta, 4 for 4th-order), and starts exploring the astrophysics of relativistic orbits and relativistic collisions.
* `inspiral.{tex,pdf}`: The fourth problem set that introduces energy decay and relativistic inspiral, and starts to explore Newtonian conditions for hydrostatic equilibrium.
* `white_dwarf.{tex,pdf}`: The fifth problem set that introduces relativistic effects in white dwarf interiors, and finding the Chandrasekhar mass limit.
* `scripts/`: A directory containing Python scripts related to these notes, e.g. a script that computes orbital separation and tidal forces as a function of gravitational wave frequency.
* `scripts/matplotlibrc`: A local configuration of matplotlib settings, because I am very picky about how my plots look.
* `Makefile`: Running `make` will produce a PDF rendering of the notes documents referred to above.

Note: all the details of this problem are captured in `*.pdf` documents in this directory. In sequntial order, they go:

1. `intro.pdf`
2. `keplerian_orbit.pdf`
3. `relativistic_orbit.pdf`
4. `inspiral.pdf`
5. `white_dwarf.pdf`

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
pdflatex -draftmode intro && \
pdflatex -draftmode intro && \
pdflatex intro
```

and the same for notes on stable orbits. Alternatively, you can use the Makefile, which `make`s life much simpler:

```
make intro
make keplerian_orbits
make relativistic_orbits
make inspiral
```

To re-generate drafts of all notes:

```
make
```

In the words of Captain Picard, `make` it so!
