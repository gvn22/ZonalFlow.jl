---
title: 'ZonalFlow: a spectral solver for dynamics and statistics of $\beta$-plane turbulence'
tags:
  - Julia
  - geophysical fluid dynamics
  - statistical state dynamics
  - quasilinear theory
  - cumulant expansions
  - generalized quasilinear theory
  - generalized cumulant expansions
authors:
  - name: Girish Nivarti^[first author,corresponding author] # note this makes a footnote saying 'co-first author'
    orcid: 0000-0003-0872-7098
    affiliation: 1 # (Multiple affiliations must be quoted)
  - name: J. Brad Marston^[second author] # note this makes a footnote saying 'co-first author'
    affiliation: 2
  - name: Steven M. Tobias^[third author]
    affiliation: 1
affiliations:
 - name: Department of Applied Mathematics, University of Leeds, Leeds, UK
   index: 1
 - name: Department of Physics, Brown University, Providence, USA
   index: 2
date: 15 August 2021
bibliography: paper.bib

# Optional fields if submitting to a AAS journal too, see this blog post:
# https://blog.joss.theoj.org/2018/12/a-new-collaboration-with-aas-publishing
aas-doi: 10.3847/xxxxx <- update this with the DOI from AAS once you know it.
aas-journal: Astrophysical Journal <- The name of the AAS journal.
---

# Summary

Fluid dynamics of the turbulent atmospheres of rotating planets and stars is often simplified to dynamics on the $\beta$-plane, which is a $2$D plane that includes effect of rotation and its gradient at a specified latitude. An important phenomenon of interest in such flows is the occurence of zonal jets, which are mean flows along the east-west direction. Prediction of the genesis and structure of turbulent zonal jets remains an open problem. To mitigate the computational requirements necessitated by direct numerical simulation, which solves the dynamical equations numerically, alternative direct statistical simulations approaches are being developed. The ZonalFlow package allows for DNS and DSS using cumulant expansions.

# Related work

# Statement of need

<!-- `Gala` is an Astropy-affiliated Python package for galactic dynamics. Python
enables wrapping low-level languages (e.g., C) for speed without losing
flexibility or ease-of-use in the user-interface. The API for `Gala` was
designed to provide a class-based and user-friendly interface to fast (C or
Cython-optimized) implementations of common operations such as gravitational
potential and force evaluation, orbit integration, dynamical transformations,
and chaos indicators for nonlinear dynamics. `Gala` also relies heavily on and
interfaces well with the implementations of physical units and astronomical
coordinate systems in the `Astropy` package [@astropy] (`astropy.units` and
`astropy.coordinates`).

`Gala` was designed to be used by both astronomical researchers and by
students in courses on gravitational dynamics or astronomy. It has already been
used in a number of scientific publications [@Pearson:2017] and has also been
used in graduate courses on Galactic dynamics to, e.g., provide interactive
visualizations of textbook material [@Binney:2008]. The combination of speed,
design, and support for Astropy functionality in `Gala` will enable exciting
scientific explorations of forthcoming data releases from the *Gaia* mission
[@gaia] by students and experts alike. -->

# Functionality

<!-- Single dollars ($) are required for inline mathematics e.g. $f(x) = e^{\pi/x}$

Double dollars make self-standing equations:

$$\Theta(x) = \left\{\begin{array}{l}
0\textrm{ if } x < 0\cr
1\textrm{ else}
\end{array}\right.$$

You can also use plain \LaTeX for equations
\begin{equation}\label{eq:fourier}
\hat f(\omega) = \int_{-\infty}^{\infty} f(x) e^{i\omega x} dx
\end{equation}
and refer to \autoref{eq:fourier} from text. -->

# Citations

<!-- Citations to entries in paper.bib should be in
[rMarkdown](http://rmarkdown.rstudio.com/authoring_bibliographies_and_citations.html)
format.

If you want to cite a software repository URL (e.g. something on GitHub without a preferred
citation) then you can do it with the example BibTeX entry below for @fidgit.

For a quick reference, the following citation commands can be used:
- `@author:2001`  ->  "Author et al. (2001)"
- `[@author:2001]` -> "(Author et al., 2001)"
- `[@author1:2001; @author2:2001]` -> "(Author1 et al., 2001; Author2 et al., 2002)"
 -->
# Figures

<!-- Figures can be included like this:
![Caption for example figure.\label{fig:example}](figure.png)
and referenced from text using \autoref{fig:example}.

Figure sizes can be customized by adding an optional second parameter:
![Caption for example figure.](figure.png){ width=20% } -->

# Acknowledgements

We would like to acknowledge support of funding from the European Research Council (ERC) under the European Unions Horizon2020 research and innovation programme (grant agreement no. D5S-DLV-786780).

# References
