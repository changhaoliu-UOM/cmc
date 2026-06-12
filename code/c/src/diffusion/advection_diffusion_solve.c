/* advection_diffusion_solve.c */
#include <errno.h>
#include <stdlib.h>
#include <string.h>
#include "double_array.h"

#include "cmc_error_message.h"
#include "cmc_diffusion_discrete_primal_weak_a_advective.h"

#include "diffusion_transient_discrete_mixed_weak.h"
#include "diffusion_transient_discrete_mixed_weak_solve_trapezoidal_next.h"
#include "diffusion_transient_discrete_mixed_weak_trapezoidal_loop_data.h"

#include "diffusion_transient_discrete_primal_weak.h"
#include "diffusion_transient_discrete_primal_weak_solve_trapezoidal_next.h"
#include "diffusion_transient_discrete_primal_weak_trapezoidal_loop_data.h"

#include "double_array.h"
#include "int.h"
#include "matrix_sparse.h"
#include "mesh.h"
#include "mesh_qc.h"

double * advection_diffusion_solve(
  struct mesh * m,
  struct matrix_sparse * m_cbd_dm1,
  double * m_inner_dm1,
  double * m_inner_d,
  struct diffusion_transient_discrete_mixed_weak * data_flow,
  double * m_inner_0,
  double * m_inner_1,
  struct diffusion_transient_discrete_primal_weak * data_diffusion,
  double time_step,
  int number_of_steps)
{
  int i, d, m_cn_dm1, m_cn_d, m_cn_0, m_cn_dm1_bar;
  int status_local = 0;

  double * flow_rate = NULL;
  double * dual_potential = NULL;
  double * concentration = NULL;
  double * rhs_final = NULL;
  double * flow_rate_reduced = NULL;
  double * y = NULL;

  struct diffusion_transient_discrete_mixed_weak_trapezoidal_loop_data * 
    input = NULL;
  struct diffusion_transient_discrete_primal_weak_trapezoidal_loop_data *
    input_diffusion = NULL;

  struct matrix_sparse * lhs_base = NULL;
  struct matrix_sparse * rhs_base = NULL;

  struct matrix_sparse * a_advective = NULL;
  struct matrix_sparse * lhs_advective = NULL;
  struct matrix_sparse * rhs_advective = NULL;
  struct matrix_sparse * lhs_new = NULL;
  struct matrix_sparse * rhs_new = NULL;

  d = m->dim;
  m_cn_dm1 = m->cn[d - 1];
  m_cn_d = m->cn[d];
  m_cn_0 = m->cn[0];

  /*
    Allocate full flow-rate history:
    flow_rate has size (number_of_steps + 1) * cn[d - 1]
    dual_potential has size (number_of_steps + 1) * cn[d]
  */
  flow_rate = (double *) malloc(
    sizeof(double) * (number_of_steps + 1) * m_cn_dm1);
  if (flow_rate == NULL)
  {
    cmc_error_message_position_in_code(__FILE__, __LINE__);
    cmc_error_message_malloc(
      sizeof(double) * (number_of_steps + 1) * m_cn_dm1,
      "flow_rate");
    goto error;
  }

  dual_potential = (double *) malloc(
    sizeof(double) * (number_of_steps + 1) * m_cn_d);
  if (dual_potential == NULL)
  {
    cmc_error_message_position_in_code(__FILE__, __LINE__);
    cmc_error_message_malloc(
      sizeof(double) * (number_of_steps + 1) * m_cn_d,
      "dual_potential");
    goto error;
  }

  /*
    Allocate concentration history:
    concentration has size (number_of_steps + 1) * cn[0]
  */
  concentration = (double *) malloc(
    sizeof(double) * (number_of_steps + 1) * m_cn_0);
  if (concentration == NULL)
  {
    cmc_error_message_position_in_code(__FILE__, __LINE__);
    cmc_error_message_malloc(
      sizeof(double) * (number_of_steps + 1) * m_cn_0,
      "concentration");
    goto error;
  }

  memset(concentration, 0,
    sizeof(double) * (number_of_steps + 1) * m_cn_0);

  /*
    rhs_final is a temporary vector for one primal weak solve.
  */
  rhs_final = (double *) malloc(sizeof(double) * m_cn_0);
  if (rhs_final == NULL)
  {
    cmc_error_message_position_in_code(__FILE__, __LINE__);
    cmc_error_message_malloc(sizeof(double) * m_cn_0, "rhs_final");
    goto error;
  }

  /*
    Step 1:
    Solve the transient flow problem first.
    This gives flow_rate at all time levels.
  */


  /*
    Step 2:
    Initialise the primal weak diffusion trapezoidal data.
    This gives the original diffusion lhs/rhs.
  */
  input
  = diffusion_transient_discrete_mixed_weak_trapezoidal_loop_data_initialize(
    m, m_cbd_dm1, m_inner_dm1, m_inner_d, data_flow, time_step);
  if (input == NULL)
  {
    cmc_error_message_position_in_code(__FILE__, __LINE__);
    fputs("cannot initialize constant input for loop\n", stderr);
    goto input_free;
  }

  input_diffusion =
    diffusion_transient_discrete_primal_weak_trapezoidal_loop_data_initialize(
      m,
      m_inner_0,
      m_inner_1,
      data_diffusion,
      time_step);

  if (input_diffusion == NULL)
  {
    cmc_error_message_position_in_code(__FILE__, __LINE__);
    fputs("cannot initialise input_diffusion\n", stderr);
    goto error;
  }

  /* initialization */
  memcpy(flow_rate, input->data->initial_flow_rate, sizeof(double) * m_cn_dm1);
  memcpy(dual_potential, input->data->initial_dual_potential, sizeof(double) * m_cn_d);
 /* the initial $n$ elements of $potential$ are the initial condition */
  memcpy(concentration, input_diffusion->data->initial, sizeof(double) * m_cn_0);
  // fprintf(stdout,"mixed : ");
  // diffusion_transient_discrete_mixed_weak_file_print_raw(stdout,input->data);
  // fprintf(stdout,"primal : ");
  // diffusion_transient_discrete_primal_weak_file_print_raw(stdout,input_diffusion->data);

  /*
    Save pure diffusion matrices.
    Every time step we will use:

      lhs = lhs_base + dt/2  * A_advective
      rhs = rhs_base - dt/2  * A_advective
  */
  lhs_base = matrix_sparse_copy(input_diffusion->lhs);
  if (lhs_base == NULL)
  {
    cmc_error_message_position_in_code(__FILE__, __LINE__);
    fputs("cannot copy lhs_base\n", stderr);
    goto error;
  }

  rhs_base = matrix_sparse_copy(input_diffusion->rhs);
  if (rhs_base == NULL)
  {
    cmc_error_message_position_in_code(__FILE__, __LINE__);
    fputs("cannot copy rhs_base\n", stderr);
    goto error;
  }

  /*
    Step 3:
    Time loop for advection-diffusion.
  */
 
  m_cn_dm1_bar = m_cn_dm1 - data_flow->boundary_neumann_dm1->a0;

  y = (double *) malloc(sizeof(double) * m_cn_d);
  if (y == NULL)
  {
    cmc_error_message_position_in_code(__FILE__, __LINE__);
    cmc_error_message_malloc(sizeof(double) * m_cn_d, "y");
    goto error;
  }
  flow_rate_reduced = (double *) malloc(sizeof(double) * m_cn_dm1_bar);
  if (flow_rate_reduced == NULL)
  {
    cmc_error_message_position_in_code(__FILE__, __LINE__);
    cmc_error_message_malloc(sizeof(double) * m_cn_dm1_bar,
                              "flow_rate_reduced");
    goto y_free;
  }

  for (i = 0; i < number_of_steps; ++i)
  {
    a_advective = NULL;
    lhs_advective = NULL;
    rhs_advective = NULL;
    lhs_new = NULL;
    rhs_new = NULL;
    status_local = 0;

    diffusion_transient_discrete_mixed_weak_solve_trapezoidal_next(
      flow_rate + m_cn_dm1 * (i + 1),
      dual_potential + m_cn_d * (i + 1),
      y,
      flow_rate_reduced,
      flow_rate + m_cn_dm1 * i,
      dual_potential + m_cn_d * i,
      input);

     fprintf(stdout,"flow_rate %d : ",i);
     double_array_file_print(stdout,m_cn_dm1,flow_rate + m_cn_dm1 * i,"--raw");
     fprintf(stdout,"\n");

    // fprintf(stdout,"flow_rate_reduced %d : ",i);
    // double_array_file_print(stdout,m_cn_dm1_bar,flow_rate_reduced,"--raw");
    // fprintf(stdout,"\n");

    /*
      Construct A_advective using flow_rate at time level i + 1.
    */
    cmc_diffusion_discrete_primal_weak_a_advective(
      &a_advective,
      &status_local,
      m,
      m_cbd_dm1,
      data_diffusion->pi_0,
      flow_rate + m_cn_dm1 * (i + 1));

    if (status_local)
    {
      cmc_error_message_position_in_code(__FILE__, __LINE__);
      cmc_error_message_cannot_calculate("a_advective");
      goto error;
    }

    // fprintf(stdout,"a_advective %d : \n",i);
    // matrix_sparse_file_print(stdout,a_advective,"--raw");

    /*
      lhs_advective =  dt/2 * A_advective
      rhs_advective = -dt/2 * A_advective
    */
    lhs_advective = a_advective;
    rhs_advective = matrix_sparse_copy(a_advective);
    matrix_sparse_scalar_multiply(
      a_advective,
      time_step * 0.5);

    if (lhs_advective == NULL)
    {
      cmc_error_message_position_in_code(__FILE__, __LINE__);
      fputs("cannot calculate lhs_advective\n", stderr);
      goto error;
    }

    matrix_sparse_scalar_multiply(
      rhs_advective,
      -time_step * 0.5);

    if (rhs_advective == NULL)
    {
      cmc_error_message_position_in_code(__FILE__, __LINE__);
      fputs("cannot calculate rhs_advective\n", stderr);
      goto error;
    }

    /*
      lhs_new = lhs_base + lhs_advective
      rhs_new = rhs_base + rhs_advective
    */
    lhs_new = matrix_sparse_linear_combination(
      lhs_advective,
      lhs_base,
      1.,
      1.);

    if (lhs_new == NULL)
    {
      cmc_error_message_position_in_code(__FILE__, __LINE__);
      fputs("cannot calculate lhs_new\n", stderr);
      goto error;
    }

     /* update Dirichlet rows of lhs_new by Dirichlet boundary conditions */
        /* apply Dirichlet boundary condition on matrix $lhs$ */
      matrix_sparse_set_identity_rows(lhs_new, input_diffusion->data->boundary_dirichlet);

    rhs_new = matrix_sparse_linear_combination(
      rhs_advective,
      rhs_base,
      1.,
      1.);

    if (rhs_new == NULL)
    {
      cmc_error_message_position_in_code(__FILE__, __LINE__);
      fputs("cannot calculate rhs_new\n", stderr);
      goto error;
    }

    /*
      Replace the matrices used by the primal weak solver.
      Free the old matrices first to avoid memory leak.
    */
    matrix_sparse_free(input_diffusion->lhs);
    matrix_sparse_free(input_diffusion->rhs);

    input_diffusion->lhs = lhs_new;
    input_diffusion->rhs = rhs_new;

    lhs_new = NULL;
    rhs_new = NULL;

    /*
      Solve concentration at time level i + 1.
    */
    diffusion_transient_discrete_primal_weak_solve_trapezoidal_next(
      concentration + m_cn_0 * (i + 1),
      rhs_final,
      concentration + m_cn_0 * i,
      input_diffusion);

    if (errno)
    {
      cmc_error_message_position_in_code(__FILE__, __LINE__);
      fputs("cannot solve concentration at next time step\n", stderr);
      goto error;
    }

    // matrix_sparse_free(a_advective);
    matrix_sparse_free(lhs_advective);
    matrix_sparse_free(rhs_advective);

    a_advective = NULL;
    lhs_advective = NULL;
    rhs_advective = NULL;
  }

  free(y);
  free(flow_rate_reduced);

  /*
    Normal cleanup.
    Do not free concentration because it is returned to main.
  */
  matrix_sparse_free(lhs_base);
  matrix_sparse_free(rhs_base);

 
  diffusion_transient_discrete_mixed_weak_trapezoidal_loop_data_free(
    input);
  diffusion_transient_discrete_primal_weak_trapezoidal_loop_data_free(
    input_diffusion);

  free(rhs_final);
  free(flow_rate);
  free(dual_potential);

  return concentration;

y_free:
  free(y);

error:
  matrix_sparse_free(a_advective);
  matrix_sparse_free(lhs_advective);
  matrix_sparse_free(rhs_advective);
  matrix_sparse_free(lhs_new);
  matrix_sparse_free(rhs_new);

  matrix_sparse_free(lhs_base);
  matrix_sparse_free(rhs_base);

  free(rhs_final);
  free(flow_rate);
  free(dual_potential);
  free(concentration);

input_free:
  if (input != NULL)
  {
    diffusion_transient_discrete_mixed_weak_trapezoidal_loop_data_free(
      input);
  }

  if (input_diffusion != NULL)
  {
    diffusion_transient_discrete_primal_weak_trapezoidal_loop_data_free(
      input_diffusion);
  }

  return NULL;
}
