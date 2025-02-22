#include <Utilities/Timer.h>
#include <Utilities/Utilities_print.h>
#include <iostream>

#include "ConvexMPCStairsLocomotion.h"
#include "GraphSearch.h"
#include "convexMPC_interface.h"

#include "Gait.h"

//оригинальные параметры MPC+WBC
#define GAIT_PERIOD 17
#define HORIZON 14

//#define GAIT_PERIOD 16
// #define GAIT_PERIOD 34 //1000 Hz

//лучшие параметры для только MPC
// #define GAIT_PERIOD 18
// #define HORIZON 5

#define STEP_HEIGHT 0.15
#define BODY_HEIGHT 0.24

// #define SHOW_MPC_SOLVE_TIME

using namespace std;

////////////////////
// Controller
////////////////////

ConvexMPCStairsLocomotion::ConvexMPCStairsLocomotion(float _dt, int _iterations_between_mpc, be2r_cmpc_unitree::ros_dynamic_paramsConfig* parameters) : iterationsBetweenMPC(_iterations_between_mpc), horizonLength(HORIZON),
                                                                                                                   //  horizonLength(10),
                                                                                                                   dt(_dt),
                                                                                                                  //  Energy(),
                                                                                                                   trotting(GAIT_PERIOD, Vec4<int>(0, GAIT_PERIOD / 2.0, GAIT_PERIOD / 2.0, 0), Vec4<int>(GAIT_PERIOD / 2.0, GAIT_PERIOD / 2.0, GAIT_PERIOD / 2.0, GAIT_PERIOD / 2.0), "Trotting"),
                                                                                                                   trotting_copy(horizonLength, Vec4<int>(0, horizonLength / 2.0, horizonLength / 2.0, 0), Vec4<int>(horizonLength / 2.0, horizonLength / 2.0, horizonLength / 2.0, horizonLength / 2.0), "Trotting_copy"),
                                                                                                                   bounding(horizonLength, Vec4<int>(5, 5, 0, 0), Vec4<int>(4, 4, 4, 4), "Bounding"),
                                                                                                                   // bounding(horizonLength, Vec4<int>(5,5,0,0),Vec4<int>(3,3,3,3),"Bounding"),
                                                                                                                   pronking(horizonLength, Vec4<int>(0, 0, 0, 0), Vec4<int>(4, 4, 4, 4), "Pronking"), jumping(horizonLength, Vec4<int>(0, 0, 0, 0), Vec4<int>(2, 2, 2, 2), "Jumping"),
                                                                                                                   // galloping(horizonLength,
                                                                                                                   // Vec4<int>(0,2,7,9),Vec4<int>(6,6,6,6),"Galloping"),
                                                                                                                   // galloping(horizonLength,
                                                                                                                   // Vec4<int>(0,2,7,9),Vec4<int>(3,3,3,3),"Galloping"),
                                                                                                                   galloping(horizonLength, Vec4<int>(0, 2, 7, 9), Vec4<int>(4, 4, 4, 4), "Galloping"), standing(horizonLength, Vec4<int>(0, 0, 0, 0), Vec4<int>(horizonLength, horizonLength, horizonLength, horizonLength), "Standing"),
                                                                                                                   // trotRunning(horizonLength, Vec4<int>(0,5,5,0),Vec4<int>(3,3,3,3),"Trot
                                                                                                                   // Running"),
                                                                                                                   trotRunning(horizonLength, Vec4<int>(0, 5, 5, 0), Vec4<int>(4, 4, 4, 4), "Trot Running"), walking(horizonLength, Vec4<int>(0, 3, 5, 8), Vec4<int>(5, 5, 5, 5), "Walking"), walking2(horizonLength, Vec4<int>(0, 5, 5, 0), Vec4<int>(7, 7, 7, 7), "Walking2"), pacing(horizonLength, Vec4<int>(5, 0, 5, 0), Vec4<int>(5, 5, 5, 5), "Pacing"), random(horizonLength, Vec4<int>(9, 13, 13, 9), 0.4, "Flying nine thirteenths trot"), random2(horizonLength, Vec4<int>(8, 16, 16, 8), 0.5, "Double Trot")
{
  _parameters = parameters;
  // discretization of the MPC model or the MPC frequency
  dtMPC = dt * iterationsBetweenMPC;
  //dt and iterationsBetweenMPC is a controller working period and connection with MPC period, could be found at FCM_State_Locomotion
  default_iterations_between_mpc = iterationsBetweenMPC;
  printf("[Convex MPC] dt: %.3f iterations: %d horizon: %d, dtMPC: %.3f\n", dt, iterationsBetweenMPC, HORIZON, dtMPC);
  
   // void setup_problem(double dt, int horizon, double mu, double f_max)
  // mu -- friction coefficient, f_max -- force limit for the MPC QP solution
  
  // setup_problem(dtMPC, horizonLength, 0.4, 1200);
  setup_problem(dtMPC, horizonLength, 0.4, 120); // original
  // setup_problem(dtMPC, horizonLength, 0.4, 650); // DH
  rpy_comp[0] = 0;
  rpy_comp[1] = 0;
  rpy_comp[2] = 0;
  rpy_int[0] = 0;
  rpy_int[1] = 0;
  rpy_int[2] = 0;

  for (int i = 0; i < 4; i++)
  {
    firstSwing[i] = true;
  }

  initSparseMPC();

  pBody_des.setZero();
  vBody_des.setZero();
  aBody_des.setZero();

  _pub_des_traj[0] = _nh.advertise<nav_msgs::Path>("/des_traj_0", 1);
  _pub_des_traj[1] = _nh.advertise<nav_msgs::Path>("/des_traj_1", 1);
  _pub_des_traj[2] = _nh.advertise<nav_msgs::Path>("/des_traj_2", 1);
  _pub_des_traj[3] = _nh.advertise<nav_msgs::Path>("/des_traj_3", 1);

  _vis_pub[0] = _nh.advertise<visualization_msgs::Marker>("/visualization_marker_0", 1);
  _vis_pub[1] = _nh.advertise<visualization_msgs::Marker>("/visualization_marker_1", 1);
  _vis_pub[2] = _nh.advertise<visualization_msgs::Marker>("/visualization_marker_2", 1);
  _vis_pub[3] = _nh.advertise<visualization_msgs::Marker>("/visualization_marker_3", 1);
}

void ConvexMPCStairsLocomotion::initialize()
{
  for (int i = 0; i < 4; i++)
  {
    firstSwing[i] = true;
  }

  firstRun = true;
}

void ConvexMPCStairsLocomotion::recompute_timing(int iterations_per_mpc)
{
  iterationsBetweenMPC = iterations_per_mpc;
  dtMPC = dt * iterations_per_mpc;
}

void ConvexMPCStairsLocomotion::_SetupCommand(ControlFSMData<float>& data)
{

  _body_height = BODY_HEIGHT;

  float x_vel_cmd, y_vel_cmd;
  float filter(0.1);

  _yaw_turn_rate = data._desiredStateCommand->rightAnalogStick[0];
  x_vel_cmd = data._desiredStateCommand->leftAnalogStick[1];
  y_vel_cmd = data._desiredStateCommand->leftAnalogStick[0];
  // _yaw_turn_rate = 0;
  // x_vel_cmd = 0;
  // y_vel_cmd = 0;

  _x_vel_des = _x_vel_des * (1 - filter) + x_vel_cmd * filter;
  _y_vel_des = _y_vel_des * (1 - filter) + y_vel_cmd * filter;
  _z_vel_des = 0.05/dtMPC*sqrt(_x_vel_des*_x_vel_des + _y_vel_des*_y_vel_des);
  _yaw_des = data._stateEstimator->getResult().rpy[2] + dt * _yaw_turn_rate;
  _roll_des = 0.;
  _pitch_des = 0.;

  // Update PD coefs
  // for sim
  //Kp = _parameters->Swing_Kp_cartesian.cast<float>().asDiagonal();
  //Kp_stance = Kp;

  //Kd = _parameters->Swing_Kd_cartesian.cast<float>().asDiagonal();
  //Kd_stance = Kd;

  // for real
  //  Kp << 150, 0, 0, 0, 150, 0, 0, 0, 150;
  //  Kp_stance = 0 * Kp;

  //  Kd << 3, 0, 0, 0, 3, 0, 0, 0, 3;
  //  Kd_stance = Kd;
}

template <>
void ConvexMPCStairsLocomotion::run(ControlFSMData<float>& data)
{
  bool omniMode = false;

  // Command Setup
  _SetupCommand(data);
  gaitNumber = data.userParameters->cmpc_gait;
  if (gaitNumber >= 10)
  {
    gaitNumber -= 10;
    omniMode = true;
  }

  auto& seResult = data._stateEstimator->getResult();

  // Check if transition to standing
  if (((gaitNumber == 4) && current_gait != 4) || firstRun)
  {
    stand_traj[0] = seResult.position[0];
    stand_traj[1] = seResult.position[1];
    stand_traj[2] = 0.21;
    stand_traj[3] = 0;
    stand_traj[4] = 0;
    stand_traj[5] = seResult.rpy[2];
    world_position_desired[0] = stand_traj[0];
    world_position_desired[1] = stand_traj[1];
  }

  // cout << "[ConvexMPCStairsLocomotion] get res done" << endl;

  // pick gait
  Gait* gait = &trotting;
  if (gaitNumber == 1)
  {
    gait = &bounding;
  }
  else if (gaitNumber == 2)
  {
    gait = &pronking;
  }
  else if (gaitNumber == 3)
  {
    gait = &random;
  }
  else if (gaitNumber == 4)
  {
    gait = &standing;
  }
  else if (gaitNumber == 5)
  {
    gait = &trotRunning;
  }
  else if (gaitNumber == 6)
  {
    gait = &random2;
  }
  else if (gaitNumber == 7)
  {
    gait = &random2;
  }
  else if (gaitNumber == 8)
  {
    gait = &pacing;
  }
  current_gait = gaitNumber;
  *gait = trotting_copy;

  //  gait->restoreDefaults();
  gait->setIterations(iterationsBetweenMPC, iterationCounter);
  //  gait->earlyContactHandle(seResult.contactSensor, iterationsBetweenMPC, iterationCounter);
  //  std::cout << "iterationCounter " << iterationCounter << std::endl;

  jumping.setIterations(iterationsBetweenMPC, iterationCounter);

  jumping.setIterations(27 / 2, iterationCounter);

  // printf("[%d] [%d]\n", jumping.get_current_gait_phase(),
  // gait->get_current_gait_phase());
  // check jump trigger
  jump_state.trigger_pressed(jump_state.should_jump(jumping.getCurrentGaitPhase()),
                             data._desiredStateCommand->trigger_pressed);

  // bool too_high = seResult.position[2] > 0.29;
  // check jump action
  if (jump_state.should_jump(jumping.getCurrentGaitPhase()))
  {
    gait = &jumping;
    recompute_timing(27 / 2);
    _body_height = _body_height_jumping;
    currently_jumping = true;
  }
  else
  {
    recompute_timing(default_iterations_between_mpc);
    currently_jumping = false;
  }

  if (_body_height < 0.02)
  {
    _body_height = BODY_HEIGHT;
    // _body_height = 0.29; //original
  }

  // integrate position setpoint
  
  Vec3<float> v_des_robot(_x_vel_des, _y_vel_des, _z_vel_des);
  Vec3<float> v_des_world = omniMode ? v_des_robot : seResult.rBody.transpose() * v_des_robot;
  Vec3<float> v_robot = seResult.vWorld;

  // pretty_print(v_des_world, std::cout, "v des world");

  // Integral-esque pitche and roll compensation
  if (fabs(v_robot[0]) > .2) // avoid dividing by zero
  {
    rpy_int[1] += dt * (_pitch_des - seResult.rpy[1]) / v_robot[0];
  }
  if (fabs(v_robot[1]) > 0.1)
  {
    rpy_int[0] += dt * (_roll_des - seResult.rpy[0]) / v_robot[1];
  }

  // cout << "cmpc se r: " << seResult.rpy[0] << " se p: " << seResult.rpy[1] <<
  // " se y: " << seResult.rpy[2] << endl;

  rpy_int[0] = fminf(fmaxf(rpy_int[0], -.25), .25);
  rpy_int[1] = fminf(fmaxf(rpy_int[1], -.25), .25);
  rpy_comp[1] = v_robot[0] * rpy_int[1];
  rpy_comp[0] = v_robot[1] * rpy_int[0] * (gaitNumber != 8); // turn off for pronking

  for (int i = 0; i < 4; i++)
  {
    pFoot[i] = seResult.position + seResult.rBody.transpose() * (data._quadruped->getHipLocation(i) + data._legController->datas[i].p);
  }

  if (gait != &standing)
  {
    world_position_desired += dt * Vec3<float>(v_des_world[0], v_des_world[1], 0);
  }

  // some first time initialization
  if (firstRun)
  {
    world_position_desired[0] = seResult.position[0];
    world_position_desired[1] = seResult.position[1];
    world_position_desired[2] = seResult.rpy[2];

    for (int i = 0; i < 4; i++)
    {
      footSwingTrajectories[i].setHeight(Vec3<float>(0,0.05,_parameters->Swing_traj_height));
      footSwingTrajectories[i].setInitialPosition(pFoot[i]);
      footSwingTrajectories[i].setFinalPosition(pFoot[i]);
    }

    firstRun = false;
  }

  // foot placement
  for (int l = 0; l < 4; l++)
  {
    swingTimes[l] = gait->getCurrentSwingTime(dtMPC, l);
  }

  float side_sign[4] = {-1, 1, -1, 1};
  float interleave_y[4] = {-0.08, 0.08, 0.02, -0.02};
  // float interleave_gain = -0.13;
  float interleave_gain = -0.2;
  // float v_abs = std::fabs(seResult.vBody[0]);
  float v_abs = std::fabs(v_des_robot[0]);

  for (int i = 0; i < 4; i++)
  {
    if (firstSwing[i])
    {
      swingTimeRemaining[i] = swingTimes[i];
    }
    else
    {
      swingTimeRemaining[i] -= dt;
    }

    footSwingTrajectories[i].setHeight(Vec3<float>(0,0.05,_parameters->Swing_traj_height));

    Vec3<float> offset(0, side_sign[i] * data._quadruped->_abadLinkLength, 0);

    Vec3<float> pRobotFrame = (data._quadruped->getHipLocation(i) + offset);

    pRobotFrame[1] += interleave_y[i] * v_abs * interleave_gain;
    float stance_time = gait->getCurrentStanceTime(dtMPC, i);

    Vec3<float> pYawCorrected = coordinateRotation(CoordinateAxis::Z, -_yaw_turn_rate * stance_time / 2) * pRobotFrame;

    Vec3<float> des_vel;
    des_vel[0] = _x_vel_des;
    des_vel[1] = _y_vel_des;
    des_vel[2] = 0.0;

    Vec3<float> Pf = seResult.position + seResult.rBody.transpose() * (pYawCorrected + des_vel * swingTimeRemaining[i]);

    //+ seResult.vWorld * swingTimeRemaining[i];

    // float p_rel_max = 0.35f;
    float p_rel_max = 0.3f;

    // Using the estimated velocity is correct
    // Vec3<float> des_vel_world = seResult.rBody.transpose() * des_vel;
    float pfx_rel = seResult.vWorld[0] * (.5 + _parameters->cmpc_bonus_swing) * stance_time * dtMPC +
                    .03f * (seResult.vWorld[0] - v_des_world[0]) +
                    (0.5f * seResult.position[2] / 9.81f) * (seResult.vWorld[1] * _yaw_turn_rate);

    float pfy_rel = seResult.vWorld[1] * .5 * stance_time * dtMPC +
                    .03f * (seResult.vWorld[1] - v_des_world[1]) +
                    (0.5f * seResult.position[2] / 9.81f) * (-seResult.vWorld[0] * _yaw_turn_rate);
    pfx_rel = fminf(fmaxf(pfx_rel, -p_rel_max), p_rel_max);
    pfy_rel = fminf(fmaxf(pfy_rel, -p_rel_max), p_rel_max);
    Pf[0] += pfx_rel;
    Pf[1] += pfy_rel;
    Pf[2] = 0.0;
   // cout << "Foot " << i << " final position desired: " << Pf.transpose() << "\n";
    footSwingTrajectories[i].setFinalPosition(Pf);
  }

  
  // calc gait
  iterationCounter++;

  // gait
  Vec4<float> contactStates = gait->getContactState();
  Vec4<float> swingStates = gait->getSwingState();
  int* mpcTable = gait->getMpcTable();

  updateMPCIfNeeded(mpcTable, data, omniMode);

  //  StateEstimator* se = hw_i->state_estimator;
  Vec4<float> se_contactState(0, 0, 0, 0);

  // ROS_INFO_STREAM("is contact: " << se_contactState(0));

  // static bool is_stance[4] = {0, 0, 0, 0};

  static nav_msgs::Path path[4];
  static geometry_msgs::PoseStamped pose[4];

  for (int foot = 0; foot < 4; foot++)
  {
    float contactState = contactStates[foot];
    float swingState = swingStates[foot];

    // if ((se_contactState(foot) == 1) && (swingState > 0) && (is_stance[foot]
    // == 0)) if ((se_contactState(foot) == 2) && (swingState > 0))
    // {
    //   swingState = 1;
    //   is_stance[foot] = 2;
    //   ROS_INFO_STREAM("Foot " << foot << " in contact early: " <<
    //   swingState);
    // }

    // if(foot == 1)
    // {
    //   ROS_INFO_STREAM("contact: " << contactState);
    //   ROS_INFO_STREAM("swing: " << swingState);
    // }

    // contactState = data._stateEstimator->getResult().contactEstimate[foot];
    // swingState = 1 - data._stateEstimator->getResult().contactEstimate[foot];

    if (swingState > 0) // foot is in swing
    {
      if (firstSwing[foot])
      {
        firstSwing[foot] = false;
        footSwingTrajectories[foot].setInitialPosition(pFoot[foot]);
        // is_stance[foot] = 0;

        path[foot].poses.clear();
        geometry_msgs::PoseStamped Emptypose;
        pose[foot] = Emptypose;
      }

      footSwingTrajectories[foot].computeSwingTrajectoryModified(swingState, swingTimes[foot],1);

      //      footSwingTrajectories[foot]->updateFF(hw_i->leg_controller->leg_datas[foot].q,
      //                                          hw_i->leg_controller->leg_datas[foot].qd,
      //                                          0); // velocity dependent
      //                                          friction compensation todo
      //                                          removed
      // hw_i->leg_controller->leg_datas[foot].qd,
      // fsm->main_control_settings.variable[2]);

      Vec3<float> pDesFootWorld = footSwingTrajectories[foot].getPosition();
      Vec3<float> vDesFootWorld = footSwingTrajectories[foot].getVelocity();
      Vec3<float> pDesLeg = seResult.rBody * (pDesFootWorld - seResult.position) -
                            data._quadruped->getHipLocation(foot);
      Vec3<float> vDesLeg = seResult.rBody * (vDesFootWorld - seResult.vWorld);
      // cout << "Foot " << foot << " relative position desired: " << pDesLeg.transpose() << "\n";
      // cout << "Foot " << foot << " relative velocity desired: " << vDesLeg.transpose() << "\n";
      //for RViz
      pose[foot].pose.position.x = pDesFootWorld.x();
      pose[foot].pose.position.y = pDesFootWorld.y();
      pose[foot].pose.position.z = pDesFootWorld.z();

      pose[foot].pose.orientation.x = 0;
      pose[foot].pose.orientation.y = 0;
      pose[foot].pose.orientation.z = 0;
      pose[foot].pose.orientation.w = 1;

      path[foot].poses.push_back(pose[foot]);

      path[foot].header.stamp = ros::Time::now();
      path[foot].header.frame_id = "odom";

      _pub_des_traj[foot].publish(path[foot]);

      // Update for WBC
      pFoot_des[foot] = pDesFootWorld;
      vFoot_des[foot] = vDesFootWorld;
      aFoot_des[foot] = footSwingTrajectories[foot].getAcceleration();

      if (!data.userParameters->use_wbc)
      {
        // Update leg control command regardless of the usage of WBIC
        data._legController->commands[foot].pDes = pDesLeg;
        data._legController->commands[foot].vDes = vDesLeg;
        data._legController->commands[foot].kpCartesian = Kp;
        data._legController->commands[foot].kdCartesian = Kd;
      }

      std::string names[4] = {"FR_hip", "FL_hip", "RR_hip", "RL_hip"};

      marker[foot].header.frame_id = names[foot];
      marker[foot].header.stamp = ros::Time();
      marker[foot].ns = "my_namespace";
      marker[foot].id = 0;
      marker[foot].type = visualization_msgs::Marker::ARROW;
      marker[foot].action = visualization_msgs::Marker::ADD;
      // pose and orientation must be zero, except orientation.w = 1
      marker[foot].pose.position.x = 0;
      marker[foot].pose.position.y = 0;
      marker[foot].pose.position.z = 0;
      marker[foot].pose.orientation.x = 0.0;
      marker[foot].pose.orientation.y = 0.0;
      marker[foot].pose.orientation.z = 0.0;
      marker[foot].pose.orientation.w = 1.0;
      marker[foot].scale.x = 0.005; // shaft diameter
      marker[foot].scale.y = 0.01;  // head diameter
      marker[foot].scale.z = 0.0;   // if not zero, specifies head length
      marker[foot].color.a = 1.0;   // Don't forget to set the alpha!
      marker[foot].color.r = 1.0;
      marker[foot].color.g = 0.0;
      marker[foot].color.b = 0.0;
      geometry_msgs::Point p1, p2;
      // start point
      p1.x = 0;
      p1.y = 0;
      p1.z = 0;
      // finish point
      p2.x = 0;
      p2.y = 0;
      p2.z = 0;
      marker[foot].points.clear();
      marker[foot].points.push_back(p1);
      marker[foot].points.push_back(p2);

      _vis_pub[foot].publish(marker[foot]);
    }
    else // foot is in stance
    {
      firstSwing[foot] = true;

      Vec3<float> pDesFootWorld = footSwingTrajectories[foot].getPosition();
      Vec3<float> vDesFootWorld = footSwingTrajectories[foot].getVelocity();
      Vec3<float> pDesLeg = seResult.rBody * (pDesFootWorld - seResult.position) - data._quadruped->getHipLocation(foot);
      Vec3<float> vDesLeg = seResult.rBody * (vDesFootWorld - seResult.vWorld);
      //cout << "Foot " << foot << " relative position desired: " << pDesLeg.transpose() << "\n";
      //cout << "Foot " << foot << " relative velocity desired: " << vDesLeg.transpose() << "\n";

      if (!data.userParameters->use_wbc) // wbc off
      {
        data._legController->commands[foot].pDes = pDesLeg;
        data._legController->commands[foot].vDes = vDesLeg;
        data._legController->commands[foot].kpCartesian = Kp_stance;
        data._legController->commands[foot].kdCartesian = Kd_stance;

        data._legController->commands[foot].forceFeedForward = f_ff[foot];
        data._legController->commands[foot].kdJoint = Vec3<float>(_parameters->Kd_joint_0, _parameters->Kd_joint_1, _parameters->Kd_joint_2).asDiagonal();

        //      footSwingTrajectories[foot]->updateFF(hw_i->leg_controller->leg_datas[foot].q,
        //                                          hw_i->leg_controller->leg_datas[foot].qd,
        //                                          0); todo removed
        // hw_i->leg_controller->leg_commands[foot].tau_ff +=
        // 0*footSwingController[foot]->getTauFF();
      }
      else
      { // Stance foot damping
        data._legController->commands[foot].pDes = pDesLeg;
        data._legController->commands[foot].vDes = vDesLeg;
        data._legController->commands[foot].kpCartesian = 0. * Kp_stance;
        data._legController->commands[foot].kdCartesian = Kd_stance;
      }
      //            cout << "Foot " << foot << " force: " <<
      //            f_ff[foot].transpose() << "\n";
      se_contactState[foot] = contactState;

      // Update for WBC
      // Fr_des[foot] = -f_ff[foot];
      std::string names[4] = {"FR_hip", "FL_hip", "RR_hip", "RL_hip"};
      marker[foot].header.frame_id = names[foot];
      marker[foot].header.stamp = ros::Time();
      marker[foot].ns = "my_namespace";
      marker[foot].id = 0;
      marker[foot].type = visualization_msgs::Marker::ARROW;
      marker[foot].action = visualization_msgs::Marker::ADD;
      // pose and orientation must be zero, except orientation.w = 1
      marker[foot].pose.position.x = 0;
      marker[foot].pose.position.y = 0;
      marker[foot].pose.position.z = 0;
      marker[foot].pose.orientation.x = 0.0;
      marker[foot].pose.orientation.y = 0.0;
      marker[foot].pose.orientation.z = 0.0;
      marker[foot].pose.orientation.w = 1.0;
      marker[foot].scale.x = 0.005; // shaft diameter
      marker[foot].scale.y = 0.01;  // head diameter
      marker[foot].scale.z = 0.0;   // if not zero, specifies head length
      marker[foot].color.a = 0.8;   // Don't forget to set the alpha!
      marker[foot].color.r = 1.0;
      marker[foot].color.g = 0.0;
      marker[foot].color.b = 0.0;
      geometry_msgs::Point p1, p2;
      // start point
      p1.x = pDesLeg[0];
      p1.y = pDesLeg[1];
      p1.z = pDesLeg[2];
      // finish point
      float koef = 500;
      p2.x = pDesLeg[0] + (-f_ff[foot][0] / koef);
      p2.y = pDesLeg[1] + (-f_ff[foot][1] / koef);
      p2.z = pDesLeg[2] + (-f_ff[foot][2] / koef);
      marker[foot].points.clear();
      marker[foot].points.push_back(p1);
      marker[foot].points.push_back(p2);
      _vis_pub[foot].publish(marker[foot]);
    }
  }

  data._stateEstimator->setContactPhase(se_contactState);
  data._stateEstimator->setSwingPhase(gait->getSwingState());

  // Update For WBC
  pBody_des[0] = world_position_desired[0];
  pBody_des[1] = world_position_desired[1];
  pBody_des[2] = _body_height;

  vBody_des[0] = v_des_world[0];
  vBody_des[1] = v_des_world[1];
  vBody_des[2] = 0.;

  aBody_des.setZero();

  pBody_RPY_des[0] = 0.;
  pBody_RPY_des[1] = 0.;
  pBody_RPY_des[2] = _yaw_des;

  vBody_Ori_des[0] = 0.;
  vBody_Ori_des[1] = 0.;
  vBody_Ori_des[2] = _yaw_turn_rate;

  contact_state = gait->getContactState();
  // END of WBC Update
}

template <>
void ConvexMPCStairsLocomotion::run(ControlFSMData<double>& data)
{
  (void)data;
  printf("call to old CMPC with double!\n");
}

void ConvexMPCStairsLocomotion::updateMPCIfNeeded(int* mpcTable, ControlFSMData<float>& data, bool omniMode)
{
  // iterationsBetweenMPC = 30;
  if ((iterationCounter % iterationsBetweenMPC) == 0)
  {
    auto seResult = data._stateEstimator->getResult();
    float* p = seResult.position.data();

    Vec3<float> v_des_robot(_x_vel_des, _y_vel_des, 0);
    Vec3<float> v_des_world = omniMode ? v_des_robot : seResult.rBody.transpose() * v_des_robot;
    // float trajInitial[12] = {0,0,0, 0,0,.25, 0,0,0,0,0,0};

    // printf("Position error: %.3f, integral %.3f\n", pxy_err[0],
    // x_comp_integral);

    // Stand gait
    if (current_gait == 4)
    {
      float trajInitial[12] = {
          _roll_des,
          _pitch_des /*-hw_i->state_estimator->se_ground_pitch*/,
          (float)stand_traj[5] /*+(float)stateCommand->data.stateDes[11]*/,
          (float)stand_traj[0] /*+(float)fsm->main_control_settings.p_des[0]*/,
          (float)stand_traj[1] /*+(float)fsm->main_control_settings.p_des[1]*/,
          (float)_body_height /*fsm->main_control_settings.p_des[2]*/,
          0,
          0,
          0,
          0,
          0,
          0};

      for (int i = 0; i < horizonLength; i++)
        for (int j = 0; j < 12; j++)
          trajAll[12 * i + j] = trajInitial[j];
    }
    else
    {
      const float max_pos_error = .1;
      float xStart = world_position_desired[0];
      float yStart = world_position_desired[1];

      if (xStart - p[0] > max_pos_error)
        xStart = p[0] + max_pos_error;
      if (p[0] - xStart > max_pos_error)
        xStart = p[0] - max_pos_error;

      if (yStart - p[1] > max_pos_error)
        yStart = p[1] + max_pos_error;
      if (p[1] - yStart > max_pos_error)
        yStart = p[1] - max_pos_error;

      world_position_desired[0] = xStart;
      world_position_desired[1] = yStart;

      float trajInitial[12] = {(float)rpy_comp[0], // 0
                               (float)rpy_comp[1], // 1
                               _yaw_des,           // 2
                               // yawStart,    // 2
                               xStart,              // 3
                               yStart,              // 4
                               (float)_body_height, // 5
                               0,                   // 6
                               0,                   // 7
                               _yaw_turn_rate,      // 8
                               v_des_world[0],      // 9
                               v_des_world[1],      // 10
                               0};                  // 11

      for (int i = 0; i < horizonLength; i++)
      {
        for (int j = 0; j < 12; j++)
          trajAll[12 * i + j] = trajInitial[j];

        if (i == 0) // start at current position  TODO consider not doing this
        {
          // trajAll[3] = hw_i->state_estimator->se_pBody[0];
          // trajAll[4] = hw_i->state_estimator->se_pBody[1];
          trajAll[2] = seResult.rpy[2];
        }
        else
        {
          trajAll[12 * i + 3] = trajAll[12 * (i - 1) + 3] + dtMPC * v_des_world[0];
          trajAll[12 * i + 4] = trajAll[12 * (i - 1) + 4] + dtMPC * v_des_world[1];
          trajAll[12 * i + 2] = trajAll[12 * (i - 1) + 2] + dtMPC * _yaw_turn_rate;
        }
      }
    }
    Timer solveTimer;

    if (_parameters->cmpc_use_sparse > 0.5)
    {
      solveSparseMPC(mpcTable, data);
    }
    else
    {
      solveDenseMPC(mpcTable, data);
    }
    // printf("TOTAL SOLVE TIME: %.3f\n", solveTimer.getMs());
  }
}

void ConvexMPCStairsLocomotion::solveDenseMPC(int* mpcTable, ControlFSMData<float>& data)
{
  auto seResult = data._stateEstimator->getResult();

  // float Q[12] = {0.25, 0.25, 10, 2, 2, 20, 0, 0, 0.3, 0.2, 0.2, 0.2};

  // float Q[12] = {0.25, 0.25, 10, 2, 2, 50, 0, 0, 0.3, 0.2, 0.2, 0.1};
  // //original
  float Q[12] = {2.5, 2.5, 10, 50, 50, 100, 0, 0, 0.5, 0.2, 0.2, 0.1};

  // float Q[12] = {0.25, 0.25, 10, 2, 2, 40, 0, 0, 0.3, 0.2, 0.2, 0.2};
  float roll = seResult.rpy[0];
  float pitch = seResult.rpy[1];
  float yaw = seResult.rpy[2];
  float* weights = Q; //Matrix 12*12 of MPC Q weights, this is her diagonal values
  float alpha = 4e-5; // MPC R weights
  //float alpha = 4e-7; // make setting eventually: DH
  float* p = seResult.position.data();
  float* v = seResult.vWorld.data();
  float* w = seResult.omegaWorld.data();
  float* q = seResult.orientation.data();

  float r[12];
  for (int i = 0; i < 12; i++)
  {
    r[i] = pFoot[i % 4][i / 4] - seResult.position[i / 4];
  }

  // printf("current posistion: %3.f %.3f %.3f\n", p[0], p[1], p[2]);

  if (alpha > 1e-4)
  {
    std::cout << "Alpha was set too high (" << alpha << ") adjust to 1e-5\n";
    alpha = 1e-5;
  }

  Vec3<float> pxy_act(p[0], p[1], 0);
  Vec3<float> pxy_des(world_position_desired[0], world_position_desired[1], 0);
  // Vec3<float> pxy_err = pxy_act - pxy_des;
  float pz_err = p[2] - _body_height;
  Vec3<float> vxy(seResult.vWorld[0], seResult.vWorld[1], 0);

  Timer t1;
  dtMPC = dt * iterationsBetweenMPC;
  setup_problem(dtMPC, horizonLength, 0.4, 120);
  // setup_problem(dtMPC,horizonLength,0.4,650); //DH
  update_x_drag(x_comp_integral);

  if (vxy[0] > 0.3 || vxy[0] < -0.3)
  {
    // x_comp_integral += _parameters->cmpc_x_drag * pxy_err[0] * dtMPC /
    // vxy[0];
    x_comp_integral += _parameters->cmpc_x_drag * pz_err * dtMPC / vxy[0];
  }

  // printf("pz err: %.3f, pz int: %.3f\n", pz_err, x_comp_integral);

  update_solver_settings(_parameters->jcqp_max_iter, _parameters->jcqp_rho, _parameters->jcqp_sigma,
                         _parameters->jcqp_alpha, _parameters->jcqp_terminate,
                         _parameters->use_jcqp);
  // t1.stopPrint("Setup MPC");
  // printf("MPC Setup time %f ms\n", t1.getMs());

  Timer t2;
  // cout << "dtMPC: " << dtMPC << "\n";
  // float* p, float* v, float* q, float* w, float* r, float roll, float pitch, float yaw, float* weights, float* state_trajectory, float alpha, int* gait
  update_problem_data_floats(p, v, q, w, r, roll, pitch, yaw, weights, trajAll, alpha, mpcTable);
  // t2.stopPrint("Run MPC");
  // printf("MPC Solve time %f ms\n", t2.getMs());

  for (int leg = 0; leg < 4; leg++)
  {
    Vec3<float> f;
    for (int axis = 0; axis < 3; axis++)
      f[axis] = get_solution(leg * 3 + axis);

    // printf("[%d] %7.3f %7.3f %7.3f\n", leg, f[0], f[1], f[2]);

    f_ff[leg] = -seResult.rBody * f;
    // Update for WBC
    Fr_des[leg] = f;
  }
}

void ConvexMPCStairsLocomotion::solveSparseMPC(int* mpcTable, ControlFSMData<float>& data)
{
  // X0, contact trajectory, state trajectory, feet, get result!
  (void)mpcTable;
  (void)data;
  auto seResult = data._stateEstimator->getResult();

  std::vector<ContactState> contactStates;
  for (int i = 0; i < horizonLength; i++)
  {
    contactStates.emplace_back(mpcTable[i * 4 + 0], mpcTable[i * 4 + 1], mpcTable[i * 4 + 2],
                               mpcTable[i * 4 + 3]);
  }

  for (int i = 0; i < horizonLength; i++)
  {
    for (u32 j = 0; j < 12; j++)
    {
      _sparseTrajectory[i][j] = trajAll[i * 12 + j];
    }
  }

  Vec12<float> feet;
  for (u32 foot = 0; foot < 4; foot++)
  {
    for (u32 axis = 0; axis < 3; axis++)
    {
      feet[foot * 3 + axis] = pFoot[foot][axis] - seResult.position[axis];
    }
  }

  _sparseCMPC.setX0(seResult.position, seResult.vWorld, seResult.orientation, seResult.omegaWorld);
  _sparseCMPC.setContactTrajectory(contactStates.data(), contactStates.size());
  _sparseCMPC.setStateTrajectory(_sparseTrajectory);
  _sparseCMPC.setFeet(feet);
  _sparseCMPC.run();

  Vec12<float> resultForce = _sparseCMPC.getResult();

  for (u32 foot = 0; foot < 4; foot++)
  {
    Vec3<float> force(resultForce[foot * 3], resultForce[foot * 3 + 1], resultForce[foot * 3 + 2]);
    // printf("[%d] %7.3f %7.3f %7.3f\n", foot, force[0], force[1], force[2]);
    f_ff[foot] = -seResult.rBody * force;
    Fr_des[foot] = force;
  }
}

void ConvexMPCStairsLocomotion::initSparseMPC()
{
  Mat3<double> baseInertia;
  baseInertia << 0.07, 0, 0, 0, 0.26, 0, 0, 0, 0.242;
  double mass = 9;
  double maxForce = 120;

  std::vector<double> dtTraj;
  for (int i = 0; i < horizonLength; i++)
  {
    dtTraj.push_back(dtMPC);
  }

  Vec12<double> weights;
  //Matrix 12*12 of MPC Q weights, this is her diagonal values
  weights << 0.25, 0.25, 10, 2, 2, 20, 0, 0, 0.3, 0.2, 0.2, 0.2;
  // weights << 0,0,0,1,1,10,0,0,0,0.2,0.2,0;

  _sparseCMPC.setRobotParameters(baseInertia, mass, maxForce);
  _sparseCMPC.setFriction(1.0);
  // _sparseCMPC.setFriction(0.4);
  _sparseCMPC.setWeights(weights, 4e-5); //send Q and R to MPC
  _sparseCMPC.setDtTrajectory(dtTraj);

  _sparseTrajectory.resize(horizonLength);
}
