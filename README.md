
This repository contains the scripts to perform Molecular Dynamics (MD) simulations of the polymer models discussed in the paper:

[Loop-extrusion and polymer phase-separation can co-exist at the single-molecule level to shape chromatin folding](https://www.biorxiv.org/content/10.1101/2021.11.02.466589v1)

To improve reproducibility, the present repository is provided with the following [DOI: 10.5281/zenodo.6726064](https://doi.org/10.5281/zenodo.6726064).

# Table of Contents

1.  [Structure of the folder and input files](#orgd8d06c7)
   1.  [Python scripts:](#org64e7fd0)
   2.  [Sub-folders:](#orgb9108a2)
2.  [Installing software](#org6d7c358)
3.  [Run simulations](#org93e1704)
4.  [Demo simulation](#orga3d6646)
5.  [Expected output of simulations](#org4b3a092)
6.  [Research projects](#org658b6d5)
7.  [Contributors](#org1af79e9)



<a id="orgd8d06c7"></a>

# Structure of the folder and input files

The current folder contains six Python scripts (.py) and three sub-folders.


<a id="org64e7fd0"></a>

## Python scripts:

The main script of MD simulations is `sbs-le.py`. It takes as input three dependencies, i.e., the scripts `saw.py` (that prepares the initial state of simulation), `sbs_class.py` (that prepares the context for SBS simulations), `le4hoomd.py` (that prepares the context for LE simulations, is taken from its own repository on <https://codeberg.org/ehsan/LE4hoomd>). The script `cont-dist-map-sbs-le.py` and `cont-dist-map-sbs-le2HDF.py` are model python scripts to produce, respectively in .txt and HDF formats, contact and distance maps from the output simulations (see below).


<a id="orgb9108a2"></a>

## Sub-folders:

-   `HCT116_chr21_30kb`: input files for LE and SBS simulations for the studied locus in HCT116 and corresponding python scripts to set model parameters. For model parameters see also the Methods sections of the paper.
-   `IMR90_chr21_30kb`: input files for LE and SBS simulations for the studied locus in IMR90 and corresponding python scripts to set model parameters. For model parameters see also the Methods sections of the paper.
-   `demo_sim`: folder containing, as an example, a very short demo simulation of the eLE model in IMR90.


<a id="org6d7c358"></a>

# Installing software

The HOOMD-blue Python package can be installed on Linux-64 or OSx-64 systems by the following command line:
`conda install -c conda-forge hoomd`
Please note that the scripts are written for HOOM version 2.9.x. The installation time requires few minutes. All details on prerequisites (e.g., software packages and libraries) can be found at the official website <http://glotzerlab.engin.umich.edu/hoomd-blue/>.


<a id="org93e1704"></a>

# Run simulations

To run one MD polymer simulation, use in the current folder the following command line:

    python3 sbs-le.py --settings $LOCUS_chr21_30kb.settings_$MODEL_$LOCUS_30kb.py --sim-id test$MODEL --hoomd cpu

-   `$LOCUS` identifies the cell line to model: set it to `IMR90` or `HCT116`.

-   `$MODEL` identifies the model to simulate: set it to `LE` or `eLE` or `SBS` or `LESBS`.

Thus, for example, to run one simulation of LE+SBS model in IMR90:

    python3 sbs-le.py --settings IMR90_chr21_30kb.settings_LESBS_IMR90_30kb.py --sim-id testLESBS --hoomd cpu

Roughly 10 minutes occur to generate a demo simulation up to 2e6 timesteps on normal PC. To have fully equilibrated conformations, our simulations run up to 2e8 timesteps, so the running time of one simulation is tens of hours. To produce an ensemble of polymer conformations, simply run multiple copies (e.g., 1000) of the above scripts.


<a id="orga3d6646"></a>

# Demo simulation

We provide as an example a demo run simulation of eLE model in IMR90 with very few timesteps (1.3e6 MD timesteps). To run the simulation, use in the current folder the following command line:

    python3 sbs-le.py --settings demo_sim.settings_eLE_IMR90_30kb.py --sim-id testeLE --hoomd cpu

Note: the simulations we performed in the paper ran for thousands of MD timesteps to have fully equilibrated conformations (up to 2e8 MD timesteps); this demo is an example.


<a id="org4b3a092"></a>

# Expected output of simulations

The main output of a simulation is a gsd file named `test$MODEL-dump.gsd`, which contains all details of system dynamics at each simulated timestep (such as x-y-z coordinates, particle types and bond types, as standard in MD simulations). $MODEL is, again, one of the models mentioned before. Hence, in the case of an eLE simulation, the output will be named: `testeLE-dump.gsd`.

The output gsd file can be easily analyzed by using standard python packages. As an example, we report the python script `cont-dist-map-sbs-le.py`, which computes distance and contact maps (in txt format) from the gsd output of the simulation. To run this python script, use the line:

    python3 cont-dist-map-sbs-le.py -i test$MODEL-dump.gsd -o heatmap --particles $P --t1 0.7 --t2 1.0 --step 0.003 --contact-thr 5.0 --dist TRUE~

As before, set `$MODEl` to `LE` or `eLE` or `SBS` or `LESBS`. Set `$P` to "750" in IMR90 and "930" in HCT116.

Similar command line holds for the python script `cont-dist-map-sbs-le2HDF.py`, which simply provides distance and contact map in HDF format.

Other outputs printed during the simulation:

-   `test$MODEL-dump-ctcf.txt` (in LE-based simulations): txt file with simulated CTCF tracks.
-   `test$MODEL-restart.gsd`: gsd file to resume MD run in case the simulation breaks.


<a id="org658b6d5"></a>

# Research projects

Simulations in the following works are performed by the code presented here:

-   [Loop-extrusion and polymer phase-separation can co-exist at the single-molecule level to shape chromatin folding](https://www.biorxiv.org/content/10.1101/2021.11.02.466589v1)


<a id="org1af79e9"></a>

# Contributors

<a href="https://github.com/ehsanirani/PhaseSeparation-LoopExtrusion-MD/graphs/contributors">
  <img src="https://contrib.rocks/image?repo=ehsanirani/PhaseSeparation-LoopExtrusion-MD" />
</a>

Made with [contrib.rocks](https://contrib.rocks).

