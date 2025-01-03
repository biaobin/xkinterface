#%% Job submission and cancaling

from my_eval import *

### alias for handling condor jobs, should be defined in ~/.zshrc
# alias "cq=condor_q"
# alias to_sub_all="ls *.submit         | awk '{print \"condor_submit \"\$0}'"
# alias to_del_all="condor_q | grep \"ID:\" | awk '{print(\"condor_rm \", \$3)}'"
# to_sub_all | bash -

### Parameter scan using Astra

### Define basic setting such as bunch charge and momenta
# Charge, pC
Q_total = 1000

# Momentum, MeV/c
P_gun = 6.3
P_booster = 17

### Define scan parameters
# FWHM, ps
var0 = [10]

# BSA, mm
var1 = [3.5]

# Gun phase, degree
var2 = np.linspace(-9, 0, 4)
var2 = [0]

# Booster phase, degree
var3 = np.linspace(-36, -16, 6)
var3 = [-28, -24, -20]

# Solenoid, A
var4 = np.linspace(360, 374, 8)
var4 = [368]

# Note: we also give gun and booster gradients here, but they are not used since they will be automatically determined by the given  momenta
combi = np.array([[v0, v1, 57.55, v2, 14.5, v3, v4] for v0 in var0 for v1 in var1 for v2 in var2 for v3 in var3 for v4 in var4])

kwargs = {}
kwargs.update(farm = True,
              submit = False,
              P_gun = P_gun,
              P_booster = P_booster)

for x in combi:
    print(x)
    obj_PhotoInjector(x, **kwargs)
    #post_PhotoInjector(x, Zstop = 5.28, **kwargs)
    pass
