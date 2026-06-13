#include "diffusion_transient_continuous.h"

/*
Transient flow test.

M = [0, 1]^3

Pressure boundary condition:
  x = 0: p = 1
  x = 1: p = 0
  y,z sides: no-flow

This gives a mild pressure gradient for checking the coupled
advection-diffusion calculation.
*/

static double pi_0(const double * x)
{
  return 1.;
}

static double kappa_1(const double * x)
{
  return 1.;
}

static double initial(const double * x)
{
  if (x[0] == 1.)
    return 1.;
  else
    return 0.;
}

static double source(const double * x)
{
  return 0.;
}

static int boundary_dirichlet(const double * x)
{
  return
    (x[0] == 0. || x[0] == 1.) &&
    (0. <= x[1] && x[1] <= 1.) &&
    (0. <= x[2] && x[2] <= 1.);
}

static double g_dirichlet(const double * x)
{
  if (x[0] == 1.)
    return 1.;
  else
    return 0.;
}

static int boundary_neumann(const double * x)
{
  return
    (
      ((x[1] == 0. || x[1] == 1.) && (0. <= x[2] && x[2] <= 1.)) ||
      ((x[2] == 0. || x[2] == 1.) && (0. <= x[1] && x[1] <= 1.))
    ) &&
    (0. <= x[0] && x[0] <= 1.);
}

static double g_neumann(const double * x)
{
  return 0.;
}

const struct diffusion_transient_continuous
diffusion_transient_advection_diffusion_test_fluid_flow =
{
  pi_0,
  kappa_1,
  initial,
  source,
  boundary_dirichlet,
  g_dirichlet,
  boundary_neumann,
  g_neumann
};