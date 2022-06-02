#ifndef DEBUG_H
#define DEBUG_H

#include <geometry_msgs/PoseWithCovarianceStamped.h>
#include <iostream>
#include <ros/ros.h>
#include <sensor_msgs/JointState.h>
#include <tf/transform_broadcaster.h>
#include <unitree_legged_msgs/AllLegsInfo.h>
#include <unitree_legged_msgs/BodyInfo.h>
#include <unitree_legged_msgs/StateError.h>

using std::cout;
using std::endl;

class Debug
{
public:
  Debug(ros::Time time_start);

  void updatePlot();
  void updateVisualization();
  void tfPublish();

  unitree_legged_msgs::AllLegsInfo all_legs_info = {};
  unitree_legged_msgs::BodyInfo body_info = {};
  float z_offset;

private:
  void _init();
  void _initPublishers();
  void _ground_truth_callback(const geometry_msgs::PoseWithCovarianceStampedConstPtr& msg);

  ros::NodeHandle _nh;
  const ros::Time _zero_time;
  ros::Time _time_start;

  ros::Publisher _pub_joint_states;
  ros::Publisher _pub_all_legs_info;
  ros::Publisher _pub_body_info;
  ros::Subscriber _sub_ground_truth;
  tf::TransformBroadcaster odom_broadcaster;
  tf::TransformBroadcaster world_odom_broadcaster;
  geometry_msgs::PoseWithCovarianceStamped _ground_trurh_pose;

private:
};

#endif // DEBUG_H
