# Implementation of Algorithms for the Measure of Inequality on Partially ordered data

This repository contains the practical work behindd my very theoretical Bachelor thesis.

## Goals 🎯

The thesis dealt with *partially ordered sets* and their implementation in R.
In particular, the implementation of Markov Chain methods for the computation of inequality measures on *linear extensions* of the poset.
After reading literature and a lot of theoretical work, the Bubley-Dier method was chosen [[1]](#1).
Refer to [the thesis](https://github.com/rmaganza/parsec-improved/blob/master/thesis/tesi_RiccardoMaganza.pdf) (in italian) for a lot more in-depth explanations.

## Tools 🛠

The work started from [the *parsec* package](https://cran.r-project.org/web/packages/parsec/index.html) which I shamelessly bumped in the *parsec_source* directory. (I'm sure my advisor won't mind).
The packaged was then extended with new-low level C functions allowing for fast generation of linear extensions and computing of inequality measures, and a couple of new R functions to interface the C level.

In particular, to implement the proposed algorithm [the *ransampl* library](https://jugit.fz-juelich.de/mlz/ransampl) was used to randomly generate numbers from a specific discrete distribution. I also included the library here with the appropriate license statements.

## LICENSE DISCLAIMER

The MIT License ONLY applies TO MY CODE in the root folder of the project. For the rest of the libraries please refer to the individual licenses in each of their folders.


## References

<a id="1">[1]</a> 
Bubley, Russ and Dyer, Martin. «Faster Random Generation of Linear Extensions». In: *Discrete Math*. 201.1-3 (1999), pp. 81–88
