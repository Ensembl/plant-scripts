# plant_tools

These are a few scripts and protocols used during the release cycle of Ensembl Plants.

These are the main steps in the cycle each time a new release is coming up:

* local_checkout.sh ; ensetup $ensembl_version

* Genomes (core dbs) from the previous release are copied from a stage server to a production server (p1,p2,p3)

* New genomes/annotations are loaded on a different production server and health-checked. If a new assembly is added, the corresponding annotation mut be loaded in the same db.

These scripts have been written following Dan Bolser's scripts at 
[https://github.com/EnsemblGenomes/personal-dbolser](https://github.com/EnsemblGenomes/personal-dbolser)
