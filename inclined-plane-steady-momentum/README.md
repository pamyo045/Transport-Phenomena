# Inclined Plane Steady Momentum Balance

## inclined_plane_steady_momentum.m
This code plots the solution to the inclined plane momentum balance of a falling film with an interactive slider to allow the user to change some of the key parameters.

### The assumptions:
  1. constant density and viscosity
  2. steady-state
  3. laminar flow (simple shear flow)
  4. fully developed flow
  5. newton's law of viscosity is applicable

### Derivation:

    d(tau)/dx=rho*g*cos(beta) reduces to: d^2(v_z)/dx^2=-(rho*g*cos(beta)/mu)
  
  where: 
  
    rho=density of fluid
    g=gravity
    beta=angle of inclination w.r.t the vertical axis
    mu=dynamic viscosity of fluid
    v_z=viscosity of fluid in the z-direction
    z=direction of flow (parallel to plane)
    x=direction orthogonal to plane
  
  After integrating twice:
  
    v_z=-(rho*g*cos(beta)/mu)*(x^2/2)+c1*x+c2
  
  Boundary conditions:
  
    d(v_z(x=0))/dx=0 
    v_z(x=delta)=0
  
  Final solution:
  
    v_z=(rho*g*delta^2*cos(beta)/(2*mu))*(1-(x/delta)^2)
