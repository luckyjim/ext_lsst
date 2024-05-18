# How used main script (and unique) of package:  script_createXML


Activate EDEN 2.1 environnement:

```console
source /cvmfs/euclid-dev.in2p3.fr/CentOS7/EDEN-2.1/bin/activate
```
make 

```console
cd EXT_PF1_dm4lsst
make purge; make; make install
```

now script is available in this directory : 

```console
InstallArea/x86_64-conda_cos6-gcc73-o2g/scripts/
```

but you must call it with run command available here

```console
build.x86_64-conda_cos6-gcc73-o2g/run
```

to give the EDEN environnement, so the global command is 

```console
build.x86_64-conda_cos6-gcc73-o2g/run EXT_LSST_Testing/EXT_PF1_dm4lsst/InstallArea/x86_64-conda_cos6-gcc73-o2g/scripts/script_createXML -p ../process_LWF1a/run_0909/butler_1696989/stage1 -t '_JMC_TEST2'
```
