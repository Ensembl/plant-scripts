# plant_tools

These are a few scripts and protocols used during the release cycle of Ensembl Plants.

These are the main steps in the cycle each time a new release is coming up:

* local_checkout.sh ; ensetup $ensembl_version

* Core dbs from the previous release are copied production server

* New assemblies/annotations are loaded on a different server and health-checked

Most of these scripts were written following Dan Bolser's scripts at 
[https://github.com/EnsemblGenomes/personal-dbolser](https://github.com/EnsemblGenomes/personal-dbolser)
