# hypercube-hmm
Hypercubic inference using Hidden Markov Models

Code for HyperHMM [1] to infer the transition probabilities on a hypercubic transition network given some observations of emitted signals from a hidden Markov model on that network. Visualisations of the inferred parameterised model and its summary dynamics are also performed. This is the expectation-maximisation cousin of HyperTraPS https://github.com/StochasticBiology/HyperTraPS (simple implementation here https://github.com/StochasticBiology/hypertraps-simple).

We also compare the behaviour of HyperHMM with HyperTraPS [2,3], Oncotrees [4], conjunctive Bayes networks (CBN) via TRONCO [5], and mutual hazard networks (MHN) [6].

Requirements
-------

The inference code uses the C++ `armadillo` library [7]. Other visualisations use R libraries `stringr` [8], `ggplot2` [9], `ggrepel` [10], `gridExtra` [11], and `igraph` [12].

Contents
=======

Examples
-------

You'll need to compile the C++ code. Install the `armadillo` library on your machine, then the command to compile from the Terminal may look something like
```
g++ hyperhmm.cpp -larmadillo -o hyperhmm.ce
```

Then a simple example could be run with

```
./hyperhmm.ce Data/simple_case1_L5.txt 5 10 simple1 1 1
```

The file `run-simple.R` attempts to demonstrate a simple implementation and plotting case using R. Alternatively, if you have run HyperHMM from the command line and just want to produce plots using the output, use `just-plot.R`.


Inference
-------

`hyperhmm.cpp` is the inference code. Compile with, for example,

`g++ hyperhmm.cpp -o hyperhmm.ce`

The code takes several command line arguments:

`./hyperhmm.ce [datafile] [number of features] [number of bootstrap resamples] [output file label] [cross-sectional data (0 or 1)] [simulate random walkers for each sample (0 or 1)]`

- `datafile` -- a datafile containing cross-sectional or longitudinal observations as bit strings (see below)
- `number of features` -- the number of features/traits/characters in the system (i.e. the length of a bit string in the datafile)
- `number of bootstrap resamples` -- 0 for a single point estimate; >0 runs the given number of resamples in addition
- `output file label` -- a string to label output files
- `cross-sectional data` -- if 1, data are assumed to be independent snapshots and the datafile should have one column. If 0, data are assumed to be observed transitions between two states and the datafile should have two columns (ancestor and descendant state)
- `simulate random walkers for each sample` -- random walkers are simulated on the point estimate hypercube to summarise the dynamics. If 1, they are simulated on every bootstrapped hypercube too and the summary is over all resamples. This can make the code take much longer for simple systems.

For example,
`./hyperhmm.ce Data/simple_case1_L5.txt 5 100 simple1 1 1`

The code outputs wide-format datafiles storing the mean and bootstrap standard deviation of the probability of observing a given change at a given time; also the trajectories of random walkers simulated on the inferred network.

Data
------
Synthetic and published data is in `Data`. The ovarian cancer dataset is from [13]; the tuberculosis dataset is from [14].

For cross-sectional observations, data should be provided as a single-column file where each row gives an independent snapshot observation of length L, for example

`0001`<br>
`0011`<br>
`1001`

For longitudinal observations, including those derived from estimated phylogenies, data should be provided as a two-column file where each row gives an independent transition observation between two states of length L, for example

`0001 0011`<br>
`0100 0111`<br>
`1001 1101`

The initial state is assumed to be all `0`s, and the system evolves by acquiring `1`s.

Plotting
------

Plotting is done in R, using functions in `hypercube-plots.R`. These include "bubble" plots for mean feature orderings, hypercube visualisations illustrating the complete transition network, and ordering graphs showing pairwise acquisition orderings.

R wrapping
-------

The code base includes several R scripts that "wrap" external calls to HyperHMM, including preparing data and pulling output for plots. This setup is less streamlined than fully integrating HyperHMM into the R environment, but we currently don't have resource for that. The code `hyperhmm-wrap.R` provides a function to call HyperHMM from R; `hypertraps-wrap.R` does the same for HyperTraPS.

The analyses and figures in the associated manuscript are reproduced with the various `...-start.R` and `...-retrieve.R` scripts. If you want to run these, pull the `Data/` contents into the working directory first.

Specifically, `double-approach-....R` compares HyperHMM and HyperTraPS, and `other-approaches-....R` compares HyperHMM, HyperTraPS, Oncotrees, MHN, and CBN (TRONCO).

References
=====

[1] Moen, M.T. and Johnston, I.G., 2023. HyperHMM: efficient inference of evolutionary and progressive dynamics on hypercubic transition graphs. Bioinformatics, 39(1), p.btac803.

[2] Johnston, I.G. and Williams, B.P., 2016. Evolutionary inference across eukaryotes identifies specific pressures favoring mitochondrial gene retention. Cell systems, 2(2), pp.101-111.

[3] Greenbury, S.F., Barahona, M. and Johnston, I.G., 2020. HyperTraPS: inferring probabilistic patterns of trait acquisition in evolutionary and disease progression pathways. Cell systems, 10(1), pp.39-51.

[4] Szabo,A. and Pappas,L. (2022) Oncotree: Estimating Oncogenetic Trees. R package version 0.3.4.

[5] De Sano,L. et al. (2016) TRONCO: an R package for the inference of cancer progression models from heterogeneous genomic data. Bioinformatics, 32, 1911–1913.

[6] Schill,R. et al. (2020) Modelling cancer progression using mutual hazard networks. Bioinformatics, 36, 241–249.

[7] Conrad Sanderson and Ryan Curtin. Armadillo: a template-based c++ library for linear algebra. Journal of Open Source Software, 1(2):26, 2016.

[8] Hadley Wickham. stringr: Simple, Consistent Wrappers for Common String Operations, 2019. R package version 1.4.0.

[9] Hadley Wickham. ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York, 2016.

[10] Kamil Slowikowski. ggrepel: Automatically Position Non-Overlapping Text Labels with ’ggplot2’, 2021. R package version 0.9.1.

[11] Baptiste Auguie. gridExtra: Miscellaneous Functions for ”Grid” Graphics, 2017. R package version 2.3.

[12] Gabor Csardi and Tamas Nepusz. The igraph software package for complex network research. InterJournal, Complex Systems:1695, 2006.

[13] Turid Knutsen et al. The interactive online sky/m-fish & cgh database and the entrez cancer chromosomes search database: linkage of chromosomal aberrations with the genome sequence. Genes, Chromosomes and Cancer, 44(1):52–64, 2005.

[14] Nicola Casali et al. Evolution and transmission of drug-resistant tuberculosis in a russian population. Nature genetics, 46(3):279–286, 2014.
