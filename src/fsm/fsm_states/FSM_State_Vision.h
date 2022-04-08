#ifndef FSM_STATE_VISION_H
#define FSM_STATE_VISION_H

#include "FSM_State.h"
#include "lcm_msgs/localization_lcmt.hpp"
#include "lcm_msgs/obstacle_visual_t.hpp"
#include "lcm_msgs/traversability_float_t.hpp"
#include "lcm_msgs/traversability_map_t.hpp"
#include "lcm_msgs/velocity_visual_t.hpp"
#include <Utilities/Timer.h>
#include <Utilities/ros_read_param.h>
#include <controllers/VisionMPC/VisionMPCLocomotion.h>
#include <controllers/WBC_Ctrl/LocomotionCtrl/LocomotionCtrl.hpp>
#include <controllers/convexMPC/ConvexMPCLocomotion.h>
#include <fstream>
#include <lcm/lcm-cpp.hpp>
#include <lcm_msgs/heightmap333_t.hpp>
#include <lcm_msgs/heightmap_t.hpp>
#include <lcm_msgs/heightnew_t.hpp>
#include <thread>

#include <geometry_msgs/PoseWithCovarianceStamped.h>
#include <grid_map_msgs/GridMap.h>
#include <ros/ros.h>
#include <tf2/convert.h>
#include <tf2/utils.h>

template<typename T>
class WBC_Ctrl;
template<typename T>
class LocomotionCtrlData;
/**
 *
 */
template<typename T>
class FSM_State_Vision : public FSM_State<T>
{
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  FSM_State_Vision(ControlFSMData<T>* _controlFSMData);

  // Behavior to be carried out when entering a state
  void onEnter();

  // Run the normal behavior for the state
  void run();

  // Checks for any transition triggers
  FSM_StateName checkTransition();

  // Manages state specific transitions
  TransitionData<T> transition();

  // Behavior to be carried out when exiting a state
  void onExit();

  std::string map_topic;
  std::string robot_pose_topic;

private:
  //-----ROS-------
  ros::NodeHandle _nh;
  ros::Subscriber _map_sub;
  ros::Subscriber _robot_pose_sub;

  grid_map_msgs::GridMap _grid_map;
  geometry_msgs::PoseWithCovarianceStamped _robot_pose;
  // Keep track of the control iterations
  int iter = 0;
  VisionMPCLocomotion vision_MPC;
  ConvexMPCLocomotion cMPCOld;

  WBC_Ctrl<T>* _wbc_ctrl;
  LocomotionCtrlData<T>* _wbc_data;

  Vec3<T> _ini_body_pos;
  Vec3<T> _ini_body_ori_rpy;
  Vec3<T> zero_vec3;
  Vec3<T> _global_robot_loc;
  Vec3<T> _robot_rpy;

  size_t x_size = 100;
  size_t y_size = 100;
  // 0930
  size_t x333_size = 300;
  size_t y333_size = 300;
  double grid_map_size = 0.02;
  double grid_detect_size = 0.02;
  double scale = grid_detect_size / grid_map_size;
  // 0928
  std::vector<std::vector<double>> height_map_2d;
  DMat<T> _height_map333; //= DMat<double>::Zero(x333_size, y333_size);
  // 1012
  size_t xnew_size = 101;
  size_t ynew_size = 101;
  double grid_size = 0.02;

  DMat<T> _height_mapnew; //= DMat<double>::Zero(101, 101);

  DMat<T> _height_map;
  DMat<int> _idx_map;
  DMat<float> idx_map;

  void handleHeightmapnewLCM(const lcm::ReceiveBuffer* rbuf, const std::string& chan,
                             const heightnew_t* msg);
  void handleHeightmap333LCM(const lcm::ReceiveBuffer* rbuf, const std::string& chan,
                             const heightmap333_t* msg);
  void handleHeightmapLCM(const lcm::ReceiveBuffer* rbuf, const std::string& chan,
                          const heightmap_t* msg);
  void handleIndexmapLCM(const lcm::ReceiveBuffer* rbuf, const std::string& chan,
                         const traversability_map_t* msg);
  void handleIndexmapfloatLCM(const lcm::ReceiveBuffer* rbuf, const std::string& chan,
                              const traversability_float_t* msg);
  void handleLocalization(const lcm::ReceiveBuffer* rbuf, const std::string& chan,
                          const localization_lcmt* msg);
  bool _b_localization_data = false;

  void visionLCMThread();

  void _elevMapCallback(const grid_map_msgs::GridMapConstPtr& msg);
  void _robotPoseCallback(const geometry_msgs::PoseWithCovarianceStampedConstPtr& msg);

  lcm::LCM _visionLCM;
  std::thread _visionLCMThread;

  vectorAligned<Vec3<T>> _obs_list; // loc, height
  obstacle_visual_t _obs_visual_lcm;

  void _updateStateEstimator();
  void _JPosStand();
  void _UpdateObstacle();
  void _UpdateObstacle_new();
  void _UpdateObstacle_trav();
  void _LocomotionControlStep(const Vec3<T>& vel_cmd);
  void _UpdateVelCommand(Vec3<T>& vel_cmd);
  void _RCLocomotionControl();
  void _Visualization(const Vec3<T>& des_vel);
  void _print_obstacle_list();

  Vec3<T> _target_pos;
};

#endif // FSM_STATE_LOCOMOTION_H
