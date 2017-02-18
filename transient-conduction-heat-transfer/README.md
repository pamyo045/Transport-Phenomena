# Transient Conduction Heat Transfer

## transient_conduction_heat_transfer.m
Solves the 1D unsteady-state temperature conduction heat equation in dimentionless form.

### Derivation:
  1-D transient heat conduction equation for a uniform rectangular plane:
  
    d(T)/d(t)=alpha*d^2(T)/(d(x))^2 

  Where:

    T     =temperature as a function of x and t (C)
    t     =time (s)
    x     =position (m)
    alpha =thermal diffusivity (m^2/s)
    T0    =initial temperature of the slab
    T1    =temperature imposed at the slab surfaces for time > 0
    b     =plane slab thickness divided by 2 (i.e. slab thickness = 2b, x=0 is the center)

  Dimensionless variable transformation:
  
    theta=(T1-T)/(T1-T0)  *dimensionless temperature*
    eta=x/b               *dimensionless position*
    tau=alpha*t/b^2       *dimensionless time* 

  Thus:
  
    d(theta)/d(tau)=d^2(theta)/d((eta))^2 

  Initial & boundary conditions:
  
    IC:   theta=1 when tau=0
    BC1:  theta=0 when eta=-1 for tau>0
    BC2:  theta=0 when eta=+1 for tau>0
