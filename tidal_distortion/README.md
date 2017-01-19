Newtonian Tidal Distortion
==========================

A working directory with tools that help illustrate and build intuition for the physical process of tidal distortion in a binary neutron star system, using nothing but Keplerian laws of motion.

Contents
--------

* `tidal_notes.{tex,pdf}`: A statement of the problem and some notes on how I solved it.
* `newtonian_tidal_forces.{py,pdf}`: A python script and its output. The script computes orbital separation and tidal forces as a function of gravitational wave frequency.
* `binary_diagram.{tex,pdf}`: A LaTeX rendering of a diagram of a neutron star binary. Pretty much does what it says on the tin.
* `matplotlibrc`: A local configuration of matplotlib settings, because I am very picky about how my plots look.
* `Makefile`: Running `make` will produce a PDF rendering of `tidal_notes.pdf`, the notes document referred to above.
