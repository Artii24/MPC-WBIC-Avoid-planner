#include "SolverMPC.h"
#include "RobotState.h"
#include "common_types.h"
#include "convexMPC_interface.h"
#include <JCQP/QpProblem.h>
#include <Utilities/Timer.h>
#include <cmath>
#include <eigen3/Eigen/Dense>
#include <eigen3/unsupported/Eigen/MatrixFunctions>
#include <qpOASES/include/qpOASES.hpp>
#include <stdio.h>
#include <sys/time.h>

#define BIG_NUMBER 5e10 // big enough to act like infinity, small enough to avoid numerical weirdness.

RobotState rs;
using Eigen::Dynamic;
using std::cout;
using std::endl;

Matrix<fpt, Dynamic, 13> A_qp;
Matrix<fpt, Dynamic, Dynamic> B_qp;
Matrix<fpt, 13, 12> Bdt;
Matrix<fpt, 13, 13> Adt;
Matrix<fpt, 25, 25> ABc, expmm;
Matrix<fpt, Dynamic, Dynamic> S;
Matrix<fpt, Dynamic, 1> X_d;
Matrix<fpt, Dynamic, 1> U_b;
Matrix<fpt, Dynamic, Dynamic> fmat;

Matrix<fpt, Dynamic, Dynamic> qH;
Matrix<fpt, Dynamic, 1> qg;

Matrix<fpt, Dynamic, Dynamic> eye_12h;

Matrix<fpt,Dynamic,1> _U_qp;
Matrix<fpt,Dynamic,1> _U_qpT;

qpOASES::real_t* H_qpoases;
qpOASES::real_t* g_qpoases;
qpOASES::real_t* A_qpoases;
qpOASES::real_t* lb_qpoases;
qpOASES::real_t* ub_qpoases;
qpOASES::real_t* q_soln;

qpOASES::real_t* H_red;
qpOASES::real_t* g_red;
qpOASES::real_t* A_red;
qpOASES::real_t* lb_red;
qpOASES::real_t* ub_red;
qpOASES::real_t* q_red;
u8 real_allocated = 0;

char var_elim[2000];
char con_elim[2000];

mfp* get_q_soln() { return q_soln; }

s8 near_zero(fpt a) { return (a < 0.01 && a > -.01); }

s8 near_one(fpt a) { return near_zero(a - 1); }
void matrix_to_real(qpOASES::real_t* dst, Matrix<fpt, Dynamic, Dynamic> src, s16 rows, s16 cols)
{
  s32 a = 0;
  for (s16 r = 0; r < rows; r++)
  {
    for (s16 c = 0; c < cols; c++)
    {
      dst[a] = src(r, c);
      a++;
    }
  }
}

void c2qp(Matrix<fpt, 13, 13> Ac, Matrix<fpt, 13, 12> Bc, fpt dt, s16 horizon)
{
  ABc.setZero();
  ABc.block(0, 0, 13, 13) = Ac;
  ABc.block(0, 13, 13, 12) = Bc;
  ABc = dt * ABc;
  expmm = ABc.exp();
  Adt = expmm.block(0, 0, 13, 13);
  Bdt = expmm.block(0, 13, 13, 12);

#ifdef K_PRINT_EVERYTHING
  cout << "Adt: \n" << Adt << "\nBdt:\n" << Bdt << endl;
#endif

  if (horizon > 19)
  {
    throw std::runtime_error("horizon is too long!");
  }

  Matrix<fpt, 13, 13> powerMats[20];
  powerMats[0].setIdentity();

  for (int i = 1; i < horizon + 1; i++)
  {
    powerMats[i] = Adt * powerMats[i - 1];
  }

  for (s16 r = 0; r < horizon; r++)
  {
    A_qp.block(13 * r, 0, 13, 13) = powerMats[r + 1]; // Adt.pow(r+1);

    for (s16 c = 0; c < horizon; c++)
    {
      if (r >= c)
      {
        s16 a_num = r - c;
        B_qp.block(13 * r, 12 * c, 13, 12) = powerMats[a_num] /*Adt.pow(a_num)*/ * Bdt;
      }
    }
  }

#ifdef K_PRINT_EVERYTHING
  cout << "AQP:\n" << A_qp << "\nBQP:\n" << B_qp << endl;
#endif
}

void resize_qp_mats(s16 horizon)
{
  int mcount = 0;
  int h2 = horizon * horizon;
  _U_qp.resize(12 * horizon, Eigen::NoChange);
  _U_qpT.resize(12 * horizon, Eigen::NoChange);

  mcount += 12 * horizon * 1;

  A_qp.resize(13 * horizon, Eigen::NoChange);
  mcount += 13 * horizon * 1;

  B_qp.resize(13 * horizon, 12 * horizon);
  mcount += 13 * h2 * 12;

  S.resize(13 * horizon, 13 * horizon);
  mcount += 13 * 13 * h2;

  X_d.resize(13 * horizon, Eigen::NoChange);
  mcount += 13 * horizon;

  U_b.resize(20 * horizon, Eigen::NoChange);
  mcount += 20 * horizon;

  fmat.resize(20 * horizon, 12 * horizon);
  mcount += 20 * 12 * h2;

  qH.resize(12 * horizon, 12 * horizon);
  mcount += 12 * 12 * h2;

  qg.resize(12 * horizon, Eigen::NoChange);
  mcount += 12 * horizon;

  eye_12h.resize(12 * horizon, 12 * horizon);
  mcount += 12 * 12 * horizon;

  // printf("realloc'd %d floating point numbers.\n",mcount);
  mcount = 0;

  A_qp.setZero();
  B_qp.setZero();
  _U_qp.setZero();
  _U_qpT.setZero();
  S.setZero();
  X_d.setZero();
  U_b.setZero();
  fmat.setZero();
  qH.setZero();
  eye_12h.setIdentity();

  if (real_allocated)
  {

    free(H_qpoases);
    free(g_qpoases);
    free(A_qpoases);
    free(lb_qpoases);
    free(ub_qpoases);
    free(q_soln);
    free(H_red);
    free(g_red);
    free(A_red);
    free(lb_red);
    free(ub_red);
    free(q_red);
  }

  H_qpoases = (qpOASES::real_t*)malloc(12 * 12 * horizon * horizon * sizeof(qpOASES::real_t));
  mcount += 12 * 12 * h2;
  g_qpoases = (qpOASES::real_t*)malloc(12 * 1 * horizon * sizeof(qpOASES::real_t));
  mcount += 12 * horizon;
  A_qpoases = (qpOASES::real_t*)malloc(12 * 20 * horizon * horizon * sizeof(qpOASES::real_t));
  mcount += 12 * 20 * h2;
  lb_qpoases = (qpOASES::real_t*)malloc(20 * 1 * horizon * sizeof(qpOASES::real_t));
  mcount += 20 * horizon;
  ub_qpoases = (qpOASES::real_t*)malloc(20 * 1 * horizon * sizeof(qpOASES::real_t));
  mcount += 20 * horizon;
  q_soln = (qpOASES::real_t*)malloc(12 * horizon * sizeof(qpOASES::real_t));
  mcount += 12 * horizon;

  H_red = (qpOASES::real_t*)malloc(12 * 12 * horizon * horizon * sizeof(qpOASES::real_t));
  mcount += 12 * 12 * h2;
  g_red = (qpOASES::real_t*)malloc(12 * 1 * horizon * sizeof(qpOASES::real_t));
  mcount += 12 * horizon;
  A_red = (qpOASES::real_t*)malloc(12 * 20 * horizon * horizon * sizeof(qpOASES::real_t));
  mcount += 12 * 20 * h2;
  lb_red = (qpOASES::real_t*)malloc(20 * 1 * horizon * sizeof(qpOASES::real_t));
  mcount += 20 * horizon;
  ub_red = (qpOASES::real_t*)malloc(20 * 1 * horizon * sizeof(qpOASES::real_t));
  mcount += 20 * horizon;
  q_red = (qpOASES::real_t*)malloc(12 * horizon * sizeof(qpOASES::real_t));
  mcount += 12 * horizon;
  real_allocated = 1;

  // printf("malloc'd %d floating point numbers.\n",mcount);

#ifdef K_DEBUG
  printf("RESIZED MATRICES FOR HORIZON: %d\n", horizon);
#endif
}

inline Matrix<fpt, 3, 3> cross_mat(Matrix<fpt, 3, 3> I_inv, Matrix<fpt, 3, 1> r)
{
  Matrix<fpt, 3, 3> cm;
  cm << 0.f, -r(2), r(1), r(2), 0.f, -r(0), -r(1), r(0), 0.f;
  return I_inv * cm;
}

// continuous time state space matrices.
void ct_ss_mats(Matrix<fpt, 3, 3> I_world,
                fpt m,
                Matrix<fpt, 3, 4> r_feet,
                Matrix<fpt, 3, 3> R_yaw,
                Matrix<fpt, 13, 13>& A,
                Matrix<fpt, 13, 12>& B,
                float x_drag)
{
  A.setZero();
  A(3, 9) = 1.f;
  A(11, 9) = x_drag;
  A(4, 10) = 1.f;
  A(5, 11) = 1.f;

  A(11, 12) = 1.f;
  A.block(0, 6, 3, 3) = R_yaw.transpose();

  B.setZero();
  Matrix<fpt, 3, 3> I_inv = I_world.inverse();

  for (s16 b = 0; b < 4; b++)
  {
    B.block(6, b * 3, 3, 3) = cross_mat(I_inv, r_feet.col(b));
    B.block(9, b * 3, 3, 3) = Matrix<fpt, 3, 3>::Identity() / m;
  }
}

void quat_to_rpy(Quaternionf q, Matrix<fpt, 3, 1>& rpy)
{
  // edge case!
  fpt as = t_min(-2. * (q.x() * q.z() - q.w() * q.y()), .99999);
  rpy(0) = atan2(2.f * (q.x() * q.y() + q.w() * q.z()), sq(q.w()) + sq(q.x()) - sq(q.y()) - sq(q.z()));
  rpy(1) = asin(as);
  rpy(2) = atan2(2.f * (q.y() * q.z() + q.w() * q.x()), sq(q.w()) - sq(q.x()) - sq(q.y()) + sq(q.z()));
}

void print_problem_setup(problem_setup* setup)
{
  printf("DT: %.3f\n", setup->dt);
  printf("Mu: %.3f\n", setup->mu);
  printf("F_Max: %.3f\n", setup->f_max);
  printf("Horizon: %d\n", setup->horizon);
}

void print_update_data(update_data_t* update, s16 horizon)
{
  print_named_array("p", update->p, 1, 3);
  print_named_array("v", update->v, 1, 3);
  print_named_array("q", update->q, 1, 4);
  print_named_array("w", update->r, 3, 4);
  pnv("Yaw", update->yaw);
  print_named_array("weights", update->weights, 1, 12);
  print_named_array("trajectory", update->traj, horizon, 12);
  pnv("Alpha", update->alpha);
  print_named_array("gait", update->gait, horizon, 4);
}

Matrix<fpt, 13, 1> x_0;
Matrix<fpt, 3, 3> I_world;
Matrix<fpt, 13, 13> A_ct;
Matrix<fpt, 13, 12> B_ct_r;

void solve_mpc(update_data_t* update, problem_setup* setup)
{
  rs.set(update->p, update->v, update->q, update->w, update->r, update->roll, update->pitch, update->yaw);

#ifdef K_PRINT_EVERYTHING

  printf("-----------------\n");
  printf("   PROBLEM DATA  \n");
  printf("-----------------\n");
  print_problem_setup(setup);

  printf("-----------------\n");
  printf("    ROBOT DATA   \n");
  printf("-----------------\n");
  rs.print();
  print_update_data(update, setup->horizon);
#endif

  // roll pitch yaw
  Matrix<fpt, 3, 1> rpy;
  quat_to_rpy(rs.q, rpy);

  // initial state (13 state representation)
  x_0 << rpy(2), rpy(1), rpy(0), rs.p, rs.w, rs.v, -9.8f;
  I_world = rs.R_yaw * rs.I_body * rs.R_yaw.transpose(); // original
  ct_ss_mats(I_world, rs.m, rs.r_feet, rs.R_yaw, A_ct, B_ct_r, update->x_drag);

#ifdef K_PRINT_EVERYTHING
  cout << "Initial state: \n" << x_0 << endl;
  cout << "World Inertia: \n" << I_world << endl;
  cout << "A CT: \n" << A_ct << endl;
  cout << "B CT (simplified): \n" << B_ct_r << endl;
#endif

  // QP matrices
  c2qp(A_ct, B_ct_r, setup->dt, setup->horizon);

  // weights
  Matrix<fpt, 13, 1> full_weight;

  for (u8 i = 0; i < 12; i++)
  {
    full_weight(i) = update->weights[i];
  }

  full_weight(12) = 0.f;

  S.diagonal() = full_weight.replicate(setup->horizon, 1);

  // trajectory
  for (s16 i = 0; i < setup->horizon; i++)
  {
    for (s16 j = 0; j < 12; j++)
    {
      X_d(13 * i + j, 0) = update->traj[12 * i + j];
    }
  }

  // note - I'm not doing the shifting here.
  s16 k = 0;
  for (s16 i = 0; i < setup->horizon; i++)
  {
    for (s16 j = 0; j < 4; j++)
    {
      U_b(5 * k + 0) = BIG_NUMBER;
      U_b(5 * k + 1) = BIG_NUMBER;
      U_b(5 * k + 2) = BIG_NUMBER;
      U_b(5 * k + 3) = BIG_NUMBER;
      U_b(5 * k + 4) = update->gait[i * 4 + j] * setup->f_max;
      k++;
    }
  }

  fpt mu = 1.f / setup->mu;
  Matrix<fpt, 5, 3> f_block;

  f_block << mu, 0, 1.f, -mu, 0, 1.f, 0, mu, 1.f, 0, -mu, 1.f, 0, 0, 1.f;

  for (s16 i = 0; i < setup->horizon * 4; i++)
  {
    fmat.block(i * 5, i * 3, 5, 3) = f_block;
  }

  qH = 2 * (B_qp.transpose() * S * B_qp + update->alpha * eye_12h);
  qg = 2 * B_qp.transpose() * S * (A_qp * x_0 - X_d);

  QpProblem<double> jcqp(setup->horizon * 12, setup->horizon * 20);

  // use jcqp = 0
  if (update->use_jcqp == 1)
  {
    jcqp.A = fmat.cast<double>();
    jcqp.P = qH.cast<double>();
    jcqp.q = qg.cast<double>();
    jcqp.u = U_b.cast<double>();

    for (s16 i = 0; i < 20 * setup->horizon; i++)
    {
      jcqp.l[i] = 0.;
    }

    jcqp.settings.sigma = update->sigma;
    jcqp.settings.alpha = update->solver_alpha;
    jcqp.settings.terminate = update->terminate;
    jcqp.settings.rho = update->rho;
    jcqp.settings.maxIterations = update->max_iterations;
    jcqp.runFromDense(update->max_iterations, true, false);
  }
  else
  {
    matrix_to_real(H_qpoases, qH, setup->horizon * 12, setup->horizon * 12);
    matrix_to_real(g_qpoases, qg, setup->horizon * 12, 1);
    matrix_to_real(A_qpoases, fmat, setup->horizon * 20, setup->horizon * 12);
    matrix_to_real(ub_qpoases, U_b, setup->horizon * 20, 1);

    for (s16 i = 0; i < 20 * setup->horizon; i++)
    {
      lb_qpoases[i] = 0.0f;
    }

    s16 num_constraints = 20 * setup->horizon;
    s16 num_variables = 12 * setup->horizon;

    qpOASES::int_t nWSR = 100;

    int new_vars = num_variables;
    int new_cons = num_constraints;

    for (int i = 0; i < num_constraints; i++)
    {
      con_elim[i] = 0;
    }

    for (int i = 0; i < num_variables; i++)
    {
      var_elim[i] = 0;
    }

    for (int i = 0; i < num_constraints; i++)
    {
      if (!(near_zero(lb_qpoases[i]) && near_zero(ub_qpoases[i])))
      {
        continue;
      }

      double* c_row = &A_qpoases[i * num_variables];
      for (int j = 0; j < num_variables; j++)
      {
        if (near_one(c_row[j]))
        {
          new_vars -= 3;
          new_cons -= 5;
          int cs = (j * 5) / 3 - 3;
          var_elim[j - 2] = 1;
          var_elim[j - 1] = 1;
          var_elim[j] = 1;
          con_elim[cs] = 1;
          con_elim[cs + 1] = 1;
          con_elim[cs + 2] = 1;
          con_elim[cs + 3] = 1;
          con_elim[cs + 4] = 1;
        }
      }
    }

    if (1 == 1)
    {
      int var_ind[new_vars];
      int con_ind[new_cons];
      int vc = 0;

      for (int i = 0; i < num_variables; i++)
      {
        if (!var_elim[i])
        {
          if (!(vc < new_vars))
          {
            printf("BAD ERROR 1\n");
          }

          var_ind[vc] = i;
          vc++;
        }
      }

      vc = 0;

      for (int i = 0; i < num_constraints; i++)
      {
        if (!con_elim[i])
        {
          if (!(vc < new_cons))
          {
            printf("BAD ERROR 1\n");
          }
          con_ind[vc] = i;
          vc++;
        }
      }

      for (int i = 0; i < new_vars; i++)
      {
        int olda = var_ind[i];
        g_red[i] = g_qpoases[olda];
        for (int j = 0; j < new_vars; j++)
        {
          int oldb = var_ind[j];
          H_red[i * new_vars + j] = H_qpoases[olda * num_variables + oldb];
        }
      }

      for (int con = 0; con < new_cons; con++)
      {
        for (int st = 0; st < new_vars; st++)
        {
          float cval = A_qpoases[(num_variables * con_ind[con]) + var_ind[st]];
          A_red[con * new_vars + st] = cval;
        }
      }

      for (int i = 0; i < new_cons; i++)
      {
        int old = con_ind[i];
        ub_red[i] = ub_qpoases[old];
        lb_red[i] = lb_qpoases[old];
      }

      if (update->use_jcqp == 0)
      {
        Timer solve_timer;
        qpOASES::QProblem problem_red(new_vars, new_cons);
        qpOASES::Options op;
        op.setToMPC();
        op.printLevel = qpOASES::PL_NONE;
        problem_red.setOptions(op);

        // H g A lb ub lbA ubA nWSR
        int rval = problem_red.init(H_red, g_red, A_red, NULL, NULL, lb_red, ub_red, nWSR);
        (void)rval;
        int rval2 = problem_red.getPrimalSolution(q_red);
        if (rval2 != qpOASES::SUCCESSFUL_RETURN)
        {
          printf("failed to solve!\n");
        }

        vc = 0;
        for (int i = 0; i < num_variables; i++)
        {
          if (var_elim[i])
          {
            q_soln[i] = 0.0f;
          }
          else
          {
            q_soln[i] = q_red[vc];
            vc++;
          }
        }
      }
      else
      { // use jcqp == 2
        QpProblem<double> reducedProblem(new_vars, new_cons);

        reducedProblem.A = DenseMatrix<double>(new_cons, new_vars);
        int i = 0;
        for (int r = 0; r < new_cons; r++)
        {
          for (int c = 0; c < new_vars; c++)
          {
            reducedProblem.A(r, c) = A_red[i++];
          }
        }

        reducedProblem.P = DenseMatrix<double>(new_vars, new_vars);
        i = 0;
        for (int r = 0; r < new_vars; r++)
        {
          for (int c = 0; c < new_vars; c++)
          {
            reducedProblem.P(r, c) = H_red[i++];
          }
        }

        reducedProblem.q = Vector<double>(new_vars);
        for (int r = 0; r < new_vars; r++)
        {
          reducedProblem.q[r] = g_red[r];
        }

        reducedProblem.u = Vector<double>(new_cons);
        for (int r = 0; r < new_cons; r++)
        {
          reducedProblem.u[r] = ub_red[r];
        }

        reducedProblem.l = Vector<double>(new_cons);
        for (int r = 0; r < new_cons; r++)
        {
          reducedProblem.l[r] = lb_red[r];
        }

        reducedProblem.settings.sigma = update->sigma;
        reducedProblem.settings.alpha = update->solver_alpha;
        reducedProblem.settings.terminate = update->terminate;
        reducedProblem.settings.rho = update->rho;
        reducedProblem.settings.maxIterations = update->max_iterations;
        reducedProblem.runFromDense(update->max_iterations, true, false);

        vc = 0;
        for (int kk = 0; kk < num_variables; kk++)
        {
          if (var_elim[kk])
          {
            q_soln[kk] = 0.0f;
          }
          else
          {
            q_soln[kk] = reducedProblem.getSolution()[vc];
            vc++;
          }
        }
      }
    }
  }

  if (update->use_jcqp == 1)
  {
    for (int i = 0; i < 12 * setup->horizon; i++)
    {
      q_soln[i] = jcqp.getSolution()[i];
    }
  }
}

void crossMatrix(Eigen::MatrixXd& R, const Eigen::VectorXd& omega)
{
  R.setZero();
  R(0, 1) = -omega(2);
  R(0, 2) = omega(1);
  R(1, 0) = omega(2);
  R(1, 2) = -omega(0);
  R(2, 0) = -omega(1);
  R(2, 1) = omega(0);
}
Eigen::Matrix<float, 13, 1>  printcost(problem_setup* setup)
{

  for (int i=0;i<setup->horizon;i++)
  {
    for (int j=0;j<12;j++)
    {
      _U_qp(i*12+j) = q_soln[i*12+j];
    }

  }
  Eigen::Matrix<float, 13, 1> F;
  F.setZero();
  _U_qpT = qH*_U_qp+qg;
  F(12)= _U_qp.adjoint()*(qH*_U_qp/2+qg);
  // _U_qp= qH*_U_qp+qg;
  F.block(0,0,12,1) =0.5*_U_qpT.block(0,0,12,1);
  // return _U_qp.block(0,0,12,1);
  return F;
}
