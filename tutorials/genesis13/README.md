# How to create jobs?

## Modify `create_jobs.py` to create batch files for the job submission, one needs to define:
  - the root name of the jobs, `jobName`
  - the seed numbers to scan, `seeds`
  - the command to call to load modules, `cmd`
  - the input `python` file, `inputName`

## Modify `run_job4.py` to prepare input files for `Genesis1.3 Version 4`, usually, one has to define:
  - in which common folder to save the simulations vs seed number, `direc`
  - the input `Astra` filename, including the path from the simulation folder, `fname`
  - number of macroparticle per slice, `npart`
  - How to generage input for `Genesis1.3`

## How to generate input for `Genesis1.3`
  - use the namelist `beam`, this defines the rms parameters of the beam
  - use the namelist `importdist`, this load the hdf5 file with 6D coordinates
  - use the namelist `importbeam`, this load the user-prepared slice-wise 6D coordinates, also in hdf5 format
  - NOTE: if `importbeam` is used, the namelist `beam` should be removed

## How to call `Genesis1.3`
In `run_job4.py/run_genesis4` function, there is a variable `cmd`. If `Genesis1.3` is already loaded with the `module load` command, then `cmd` will be the name of the executable, like `genesis4`; if not, then one should give the full path to the executable.
