# Transport the beam with `THz-2nC-booster-to-undulator.py`

This script has defined quadrupole gradients from EMSY1 to undulator for 2 nC, one can try to run it to get familiar with `Ocelot`.

The input `Astra` beam can be found from [here](https://syncandshare.desy.de/index.php/s/kR7WsQPQSa8FbgB)


# Tuning of quadupole strengths

- `DoubletTuning.py`, give the strength of the first quadrupole and tune the second one in order to get a round beam downstream at a defined position.
- `TripletTuning.py`, give the strength of the first quadrupole, tune the second and third quadrupoles for a round beam transport downstream
- `QuadrupletTuning.py`, give the strength of the first and second quadrupoles, tune the third and fourth for a round beam transport downstream

# Matching of the beam