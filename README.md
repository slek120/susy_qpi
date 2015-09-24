# Setup

    make

# Execute

    ./susy_qpi

# Configure

#### Integration

Change minpts, maxpts, abserr, and relerr in susy_qpi.f90


#### Data

Change t, mu, x0, epsf, V, Uc, Uf, and omega in susy_qpi.f90

    make
    ./susy_qpi

# Output

susy_qpi.log keeps a log of configuration and errors

watch the log during execution with this command

    tail -f susy_qpi.log

YYYYMMDD\_HHSS\_susy\_qpi.dat has data in the format

| qx | qy | g(qx,qy,w) | estimated error |
|----|----|------------|-----------------|