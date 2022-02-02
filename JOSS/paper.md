---
title: 'dNami: a framework for solving systems of balance laws using explicit numerical schemes.'
tags:
  - Python
  - astronomy
  - dynamics
  - galactic dynamics
  - milky way
authors:
  - name: Nicolas Alferez^[co-first author] # note this makes a footnote saying 'co-first author'
    orcid: 0000-0003-3482-4529 
    affiliation: 1 
  - name: Emile Touber^[co-first author] ^[corresponding author] # note this makes a footnote saying 'co-first author'
    orcid: 0000-0002-7553-0351 
    affiliation: "2, 3"
  - name: Stephen D. Winn 
    orcid: 0000-0001-5972-4271 
    affiliation: "2, 3"
  - name: Yussuf Ali  
    affiliation: 2
affiliations:
 - name: Laboratoire DynFluid, Arts et Métiers ParisTech, Paris, France  
   index: 1
 - name: Okinawa Institute of Science and Technology, Okinawa, Japan 
   index: 2
 - name: Department of Mechnanical Engineering, Imperial College London, London, UK
   index: 3
date: 2 February 2022
bibliography: paper.bib

# Optional fields if submitting to a AAS journal too, see this blog post:
# https://blog.joss.theoj.org/2018/12/a-new-collaboration-with-aas-publishing
# aas-doi: 10.3847/xxxxx <- update this with the DOI from AAS once you know it.
aas-journal: Astrophysical Journal <- The name of the AAS journal.
---

# Summary

Solving \autoref{eq:gov_eq} 

\begin{equation} \label{eq:gov_eq}
	\frac{\partial}{\partial t} \textbf{q} = \textbf{RHS}(\textbf{q})
\end{equation}

# Statement of need


# Features 

Specification of governing equations in symbolic form. Automatic translation to FORTRAN. Flexiblility of finite difference scheme (stencil and order) and filter. Speed with FORTRAN code. Interaction with the python interface at runtime.   



# Citations

Citations to entries in paper.bib should be in
[rMarkdown](http://rmarkdown.rstudio.com/authoring_bibliographies_and_citations.html)
format.

If you want to cite a software repository URL (e.g. something on GitHub without a preferred
citation) then you can do it with the example BibTeX entry below for @fidgit.

For a quick reference, the following citation commands can be used:
- `@author:2001`  ->  "Author et al. (2001)"
- `[@author:2001]` -> "(Author et al., 2001)"
- `[@author1:2001; @author2:2001]` -> "(Author1 et al., 2001; Author2 et al., 2002)"

# Figures

Figures can be included like this:
![Caption for example figure.\label{fig:example}](figure.png)
and referenced from text using \autoref{fig:example}.

Figure sizes can be customized by adding an optional second parameter:
![Caption for example figure.](figure.png){ width=20% }

# Acknowledgements

We acknowledge contributions from Brigitta Sipocz, Syrtis Major, and Semyeong
Oh, and support from Kathryn Johnston during the genesis of this project.

# References
