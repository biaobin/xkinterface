cd Q-1000.00pC-D-3.50mm-E1-57.55MV_m-phi1-0.00deg-E2-15.10MV_m-phi2--36.00deg-I-360.00A
generator gen.in 2>&1 | tee gen.log

astra ast.in 2>&1 | tee ast.log
