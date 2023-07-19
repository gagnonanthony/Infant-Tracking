Infant-Tracking
===============
Complete pipeline to perform tractography from infant diffusion MRI data. Adapted from the SCIL TractoFlow Pipeline (https://github.com/scilus/tractoflow.git)

SINGULARITY
-----------
If you are running this pipeline on Linux, it is recommended you use the SCIL singularity container that contains all the relevant dependencies. 
You can get the image by running this command: 

`` sudo singularity pull scilus.sif docker://scilus/scilus:latest``

DOCKER
------
If you are on MacOS or Windows, you can use Docker to run Infant-Tracking.
Prebuilt image are available here:

https://hub.docker.com/r/scilus/scilus

DEPENDENCIES
------------
You can also run Infant-Tracking without any container, but you need these dependencies installed on your machine to make it work:

* Scilpy (https://github.com/scilus/scilpy.git)
* Mrtrix3 (https://www.mrtrix.org/)
* FSL (https://fsl.fmrib.ox.ac.uk/fsl/fslwiki)
* ANTs (http://stnava.github.io/ANTs/)

All these dependencies must be in your PATH environment variable in order to launch Infant-Tracking. 

USAGE
-----
See _USAGE_ or run `` nextflow run main.nf --help `` for more details.

REFERENCES
----------
This pipeline is adapted from the SCIL TractoFlow pipeline, see:

Theaud, G., Houde, J.-C., Boré, A., Rheault, F., Morency, F., Descoteaux, M.,
    TractoFlow: A robust, efficient and reproducible diffusion MRI pipeline
    leveraging Nextflow & Singularity, NeuroImage,
    https://doi.org/10.1016/j.neuroimage.2020.116889.

