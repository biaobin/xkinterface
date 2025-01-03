# What can it do

1. Define the input parameters for `Astra` and `Generator` using `**kwargs`, e.g.,
```python
# Load the modules
from interface import *
# Define Input module
generator = Generator1(FNAME = 'beam.ini', IPart = 50000, Species = 'electrons', Q_total = -4.,
                      Ref_Ekin = 0.0e-6, LE = 0.55e-3, dist_pz = 'i',
                      Dist_z = 'p', Lt = 21.5e-3, rt = 2e-3, Cathode = True,
                      Dist_x = 'r', sig_x = sigma_x, Dist_px = 'g', Nemit_x = 0,
                      Dist_y = 'r', sig_y = sigma_x, Dist_py = 'g', Nemit_y = 0)
generator.set(Q_total = -0.1, Ipart = 10000)
print(generator.output)
```

2. Generate input files for running `Astra`, including batch files for job submitting, e.g.,
```python
generator.write(direc+os.sep+gen_name+'.in')
astra.write(direc+os.sep+ast_name+'.in')

astra.submit(job_name, ast_name, gen_name, direc)
```
This will generate the input files saved in the folder `direc` and job submitting files saved in the current folder.

3. The above can be done using a script. By defing the parameters, one can generate input files for parameter scan, which is efficient for 2-4 variables.

4. By combining it with optimization algorithms, one can do global optimization with more variables.

# How to use it

There is a tutorial notebook [here](https://gitlab.desy.de/xiangkun.li/interface/-/blob/master/tutorials/astra/astra_demo.ipynb?ref_type=heads).

# Examples

- `my_object.py` defines an objective function and an post-processing function. By calling one of them, one can run `Astra` or analyze the outputs from `Astra`.

- `ParaScan.py` defines the script for parameter scan.

- `test_pymoo_MPI.py` defines the script for global optimization using `NSGA-II` from `pymoo`.


