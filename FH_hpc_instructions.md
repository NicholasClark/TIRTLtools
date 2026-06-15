# Instructions for installing and running TIRTLtools on the Fred Hutch HPC

## Installing python dependencies

### Open a terminal emulator

Open “PowerShell” on windows or “terminal” on Mac and type the following
to connect to the HPC:

    ssh maestro # (or ssh rhino02)

### Load the MiniForge module to use conda/mamba:

    ml MiniForge3 # load miniforge for conda/mamba

### Create a conda environment for TIRTLtools:

Note: the installation of cupy, the main python dependency, is dependent
on the GPU and the CUDA software version, which is different on
`maestro` and `rhino` servers. This means you will need a separate conda
environment for each server, installed when connected to the
corresponding server.

Note 2: `maestro` has faster GPUs (`chorus` partition), but is only
accessible through the terminal, not through Fred Hutch’s OnDemand
Rstudio instances. `rhino` has slower GPUs, but is available through
Rstudio. If running `TCRdist`, `run_pairing`, or `cluster_tcrs` through
Rstudio, you need to request a gpu and note the R version (current
default is `4.4.0_foss_2023b`). TIRTLtools will need to be installed on
this same R version. You can check available R versions from the
terminal with `module spider R`.

    conda create --name tirtltools_maestro # (or tirtltools_rhino)

### Activate the conda environment:

    conda init
    source ~/.bashrc
    conda activate tirtltools_maestro # (or tirtltools_rhino)

### Install the python dependencies:

Note: mamba is a faster drop-in replacement of conda. This may still
take a long time (10-30 minutes)

    mamba install scipy numpy cupy pandas

Now you are done installing the TIRTLtools python dependencies!

## Installing TIRTLtools in R

### Connect to the HPC:

    ssh rhino #(or maestro)

### Load the R module:

    ml R
    R

### Install R dependencies:

Note: this may take a long time. (10-30 minutes)

    if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
    if (!requireNamespace("reticulate", quietly = TRUE)) install.packages("reticulate")
    if (!requireNamespace("TIRTLtools", quietly = TRUE)) remotes::install_github("NicholasClark/TIRTLtools")

## Loading and using TIRTLtools on the HPC

### Connect to the HPC:

    ssh maestro #(or rhino02)

### Load the R module:

    ml R
    R

### Load your conda environment:

Note: if you don’t know the path of your conda environment, you can find
it by typing `conda env list` in the HPC terminal. If it complains about
the conda command, you need to run `source ~/.bashrc` and then run
`conda env list` again.

    reticulate::use_condaenv("~/.conda/envs/tirtltools_maestro") # path to your conda environment

If this doesn’t work, try something like this, also specifying the conda
binary:

    use_condaenv("~/.conda/envs/tirtltools_maestro", conda = "/app/software/Miniforge3/24.1.2-0/bin/conda")

### Load TIRTLtools:

    library(TIRTLtools)

### Test the GPU functionality:

Run TCRdist on the VDJ-db dataset from the package:

    tmp = TCRdist(vdj_db)

You should see a printout of the progress of the function.
