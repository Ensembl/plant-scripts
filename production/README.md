# Ensembl Plants production

These are scripts and protocols used during the release cycle of Ensembl Plants.

This is the usual order in which they are run:

* misc_scripts/local_checkout.sh ; ensetup $ensembl_version

* scripts in core/ 

* scripts in variation/

* scripts in funcgen/

<!--Most of these scripts were written following Dan Bolser's scripts at 
[https://github.com/EnsemblGenomes/personal-dbolser](https://github.com/EnsemblGenomes/personal-dbolser) -->

A few important notes:

[1] These scripts are not CI-tested; we changed as needed across releases

[2] Some scripts use environmental variables set up in our .bashrc files

[3] Some scripts use modules under the Tools directory. In order to use these you should run (or add to your .bashrc file)
PERL5LIB=<local_path_of_repo>:$PERL5LIB
For example: PERL5LIB=/nfs/production/panda/ensemblgenomes/development/gnaamati/plant_tools:$PERL5LIB
