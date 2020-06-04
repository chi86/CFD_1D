# CFD_1D
CFD code for 1D channel/pipe flow
fluid material properties based on

"C Irrenfried, Convective turbulent near wall heat transfer at high Prandtl numbers: A modelling approach based on Direct Numerical Simulations and experiments"

https://diglib.tugraz.at/convective-turbulent-near-wall-heat-transfer-at-high-prandtl-numbers-2020


## work in progress
- postprocessing routin not yet public available, all results are put into "data.dat"





# Compilation

## debug or verbose

```bash
make verbose
make debug
``` 


# Makefile modification

## Laminar flow

```make
CFLAGS += -DLAMINAR
``` 

## Double precision

```make
CFLAGS += -DUSE_DOUBLES=1
``` 

## Channel or pipe flow

```make
all: channel
```
or
```make
all: pipe
```

## temperatur dependent fluid properties

```make
all: pipe_var
```

## constant material fluid properties

```make
all: pipe
``` 

