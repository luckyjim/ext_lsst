# LSST stage 1 : User manual

The Software User Manual is available on the Redmine OU-EXT page [here](https://euclid.roe.ac.uk/projects/sgu/wiki/EXT-SUMLSST_stage1_SUM)

We summarize the main stages:
* Under LSST env
    * Collect star reference catalog in the text file
    * Collect raw image in text file
    * Create a parameters file
    * Launch script *create_jobs_2step.py*
* Outside LSST env 
    * Submit script *sge_arrayjob_lsst.sh* "LSST step" on SGE batch queue
    * Wait the end of all jobs : FITS product are created in local disk /sps
    * Submit script *sge_arrayjob_euclid.sh* "Euclid step" on SGE batch queue
    * Wait the end of all jobs : FITS product are ingested in EAS


# Description directory

## Scripts

* **create_arrayjob_lsst.py**

Create with parameters file 2 scripts :
* *sge_arrayjob_lsst.sh*
* *sge_arrayjob_euclid.sh*

for the scheduler of CCIN2P3 data center of type [array jobs](https://doc.cc.in2p3.fr/en/Computing/job-types/job-array.html) to processing all CCD of all visit in 2 steps. Others scripts present are used by these 2 main scripts:
* **butler4euclid.py** : Create a butler and ingest a reference star catalog
* **script_lsst_step.py** : stage 1 "LSST step" for one focal plane
* **script_euclid_step.py** : stage 1 "Euclid step" for one focal plane

## Modules

* **process_image_lsst.py** : main library
* file_conf.py : manage file parameters

All scripts are write with *process_image_lsst.py* library
