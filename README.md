# AutoGenes

**AutoGenes** generates worklists for the automation of the SpeedyGenes method
for assembly of protein variant libraries:

> SpeedyGenes: an improved gene synthesis method for the efficient production of
error-corrected, synthetic protein libraries for directed evolution.
Currin A, Swainston N, Day PJ, Kell DB. *Protein Eng Des Sel*. 2014, **27**:273-80.

To run AutoGenes on Windows, open a Command Prompt, navigate to this directory
and type:

`run.bat`

For Mac / Linux, open a Terminal window, navigate to this directory and type:

`bash run.sh`

Output will be generated in the `out` directory.

By default, these scripts will run **AutoGenes** according to the parameters in
the `run.bat` or `run.sh` file.

The scripts run a Python command as follows:

`python autogenes/run.py data/plates 2 3 out/ MAON`

where:

* `data/plates` specifies the directories in which the plate maps (as csv
files) are located.
    * Two plate files are required, named `wt.csv` containing
wildtype oligos, and `mut.csv` containing mutant oligos.
    * The existing files indicate the format and headings required.

* `2` specifies the maximum number of mutants to consider per oligo;

* `3` specifies the number of blocks to be assembled into a full length gene;

* `out` specifies the directory to which to write output;

* `MAON` is a short project name used to generate plate identifiers.