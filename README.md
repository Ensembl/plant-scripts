# plant_tools

These are a scripts and protocols used during the release cycle of Ensembl Plants.

This is the usual order in which they are run:

* local_checkout.sh ; ensetup $ensembl_version

* scripts in core/ 

* scripts in variation/

* scripts in funcgen/

Most of these scripts were written following Dan Bolser's scripts at 
[https://github.com/EnsemblGenomes/personal-dbolser](https://github.com/EnsemblGenomes/personal-dbolser)

* A few things to take into account:
Some of these scripts use modules in under the Tools directory
In order to use these modules you should run (or add to your .bashrc file)
PERL5LIB=<local_path_of_repo>:$PERL5LIB
For example: PERL5LIB=/nfs/production/panda/ensemblgenomes/development/gnaamati/plant_tools:$PERL5LIB
