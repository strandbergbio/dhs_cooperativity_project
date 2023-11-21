Using data from Vierstra et al 2020 to estimate cooperativity between DHS and TFBS sites.


# Requirements
Building the occupancy data requires [`tabix`](https://www.htslib.org/doc/tabix.html) and [`bgzip`](https://www.htslib.org/doc/bgzip.html), both of which are installed with [HSTLib](https://www.htslib.org/download/). It also requires an internet connection as `tabix` pulls the nucleotide-level data from [vierstra.org](https://www.vierstra.org/resources/dgf).

The scripts allow you to specify a root directory for storing data and otherwise assumes that the data is organized according to the directory structure on [vierstra.org](https://www.vierstra.org/resources/dgf).

