#pragma once

#include <Utilities/utilities.h>
#include <geometry_msgs/PoseWithCovarianceStamped.h>
#include <geometry_msgs/TransformStamped.h>
#include <iostream>
#include <nav_msgs/Odometry.h>
#include <nav_msgs/Path.h>
#include <ros/ros.h>
#include <sensor_msgs/Imu.h>
#include <sensor_msgs/JointState.h>
#include <tf/transform_broadcaster.h>
#include <tf2_ros/transform_broadcaster.h>
#include <tf2_ros/transform_listener.h>
#include <unitree_legged_msgs/AllLegsInfo.h>
#include <unitree_legged_msgs/BodyInfo.h>
#include <unitree_legged_msgs/StateError.h>
#include <unitree_legged_msgs/MetricInfo.h>
#include <nav_msgs/Odometry.h>
#include <sensor_msgs/Imu.h>
#include <cppTypes.h>

// #include <controllers/convexMPC/Metrics.h>
#include <visualization_msgs/Marker.h>

using std::cout;
using std::endl;

#define PUB_IMU_AND_ODOM

// Класс содержит структуры данных с информацией о состоянии робота
// По указателю на экземпляр этого класса из разных участков контроллера
// передается информация в структуры этого класса.
// В классе реализованы методы для создания и публикации сообщений в
// ROS, в том числе, для визуализации в rviz.
struct MetricData
{ 
  Vec4<float> final_body_cost; 
  Vec4<float> final_cost;
  Vec4<float> final_leg_cost;
  Vec3<float> gradient_cost[4];
  float mpc_cost;
};

class Debug
{
public:
  Debug(ros::Time time_start);

  void updatePlot();
  void updateVisualization();
  void tfOdomPublishRS_t265(ros::Time stamp);
  void tfOdomPublish(ros::Time stamp);
  void tfPublish();
  void updateMetrics();

  unitree_legged_msgs::AllLegsInfo all_legs_info = {};
  unitree_legged_msgs::BodyInfo body_info = {};
  float z_offset;
  nav_msgs::Odometry ground_truth_odom = {};
  sensor_msgs::Imu imu;
  geometry_msgs::Point last_p_stance[4] = {};
  geometry_msgs::Point last_p_local_stance[4] = {};
  geometry_msgs::Point mnk_plane = {};
  Vec3<float> hip_location[4] = {};
  Mat3<float> Rbody = {};

  nav_msgs::Path leg_traj_des[4];
  geometry_msgs::Point leg_force[4];

  bool is_map_upd_stop;
  ros::Time time_stamp_udp_get;
  double vio_z;
  MetricData metric_data = {};

private:
  void _init();
  void _initPublishers();
  Vec3<float> _getHipLocation(uint8_t leg_num);
  void _drawLastStancePoints();
  void _drawSwingFinalPoints();
  void _drawEstimatedStancePLane();
  void _drawLegsDesiredTrajectory();
  void _drawLegsForce();
  void _drawLocalBodyHeight();

  ros::NodeHandle _nh;
  const ros::Time _zero_time;
  ros::Time _time_start;

  ros::Subscriber _sub_ground_truth;

  ros::Publisher _pub_joint_states;
  ros::Publisher _pub_imu;
  ros::Publisher _pub_all_legs_info;
  ros::Publisher _pub_odom;
  ros::Publisher _pub_body_info;
  ros::Publisher _pub_metric_info;
  // tf::TransformBroadcaster odom_broadcaster;
  // tf::TransformBroadcaster world_broadcaster;
  ros::Publisher _pub_vis_last_p_stance;
  ros::Publisher _pub_vis_swing_pf;
  ros::Publisher _pub_vis_estimated_stance_plane;
  ros::Publisher _pub_vis_leg_des_traj[4];
  ros::Publisher _pub_vis_leg_force[4];
  ros::Publisher _pub_vis_local_body_height;
  tf2_ros::TransformBroadcaster odom_broadcaster;
  tf2_ros::TransformBroadcaster world_odom_broadcaster;
  tf2_ros::Buffer _tf_buffer;
  tf2_ros::TransformListener _tf_listener;
};
