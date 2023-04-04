# prefpp

This repository provides the C++ implementation of a preference-based postprocessing method proposed in the following paper:

> Ryoji Tanabe: **On the Unbounded External Archive and Population Size in Preference-based Evolutionary Multi-objective Optimization Using a Reference Point**, Proc. ACM Genetic and Evolutionary Computation Conference [(GECCO2023)](https://gecco-2023.sigevo.org), pdf, supplement

This repository uses the C++ implementation of the ND-Tree-based update method proposed in the following paper:

> Andrzej Jaszkiewicz, Thibaut Lust: **ND-Tree-Based Update: A Fast Algorithm for the Dynamic Nondominance Problem**. IEEE Trans. Evol. Comput. 22(5): 778-791 (2018)

I downloaded the code of the ND-Tree-based update method from [the supplemental website](https://sites.google.com/site/ndtreebasedupdate/). However, it seems that the link has expired.

I implemented the preference-based postprocessing method by revising ``Main.cpp``. I also the type of ``ObjectiveValues`` from integer to double by revising ``solution.h``. In addition, I added a new method ``saveTo2dVec`` to the class ``TTreeSet`` in ``ttreeeset.h``.

## Requirements

This code require a C++ compiler. I compiled this code with GNU g++ version 11.3.0 on the Ubuntu 22.04.

## Usage

After compiling the code by ``make``, the following example can be performed, where it will take about 1.5 minutes.

> ./Main -res_file_path sample/RNSGA2_mu100_DTLZ1_m2_z0.6_0.4_0th_run.csv -pp_file_path pp_res.csv -pref_point_file_path ref_point_data/m2.csv -n_obj 2 -pset_size 100 -roi_radius 0.1

- ``res_file_path`` is a file path to the results of a PBEMO algorithm. The above example file ``sample/RNSGA2_mu100_DTLZ1_m2_z0.6_0.4_0th_run.csv`` provides all $5 \times 10^4$ points found so far by R-NSGA-II with the population size 100 using the reference point $\mathbf{z}=(0.6, 0.4)$ on the bi-objective DTLZ1 problem. First, only non-dominated points are selected from the $5 \times 10^4$ points by the ND-Tree-based update method. Then, the preference-based postprocessing method is applied to the non-dominated point set.
- ``pp_file_path`` is a file path that saves the postprocessing result. In the above example, the postprocessing result is saved to ``pp_res.csv``.
- ``pref_point_file_path`` is a file path to the reference point. 
- ``n_obj`` specifies the number of objectives.
- ``pset_size`` specifies the subset size $k$.
- ``roi_radius`` determines the radius of the approximated ROI $r$.
