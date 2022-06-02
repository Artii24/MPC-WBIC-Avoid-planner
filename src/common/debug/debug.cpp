#include "debug.hpp"

Debug::Debug(ros::Time time_start)
  : _zero_time(0)
  , _time_start(time_start)
{
  z_offset = 0.f;
  _init();
  _sub_ground_truth = _nh.subscribe("ground_truth_pose", 1, &Debug::_ground_truth_callback, this,
                                    ros::TransportHints().tcpNoDelay(true));
}

void Debug::_init()
{
  _initPublishers();

  body_info.quat_act.w = 1.0;
}

void Debug::_ground_truth_callback(const geometry_msgs::PoseWithCovarianceStampedConstPtr& msg)
{
  _ground_trurh_pose = *msg;
}

void Debug::_initPublishers()
{
  _pub_joint_states = _nh.advertise<sensor_msgs::JointState>("/joint_states", 1);
  _pub_all_legs_info = _nh.advertise<unitree_legged_msgs::AllLegsInfo>("/all_legs_info", 1);
  _pub_body_info = _nh.advertise<unitree_legged_msgs::BodyInfo>("/body_info", 1);
  // _pub_parameters = _nh.advertise<unitree_legged_msgs::Parameters>("/parameters", 1);
}

void Debug::updatePlot()
{
  ros::Duration delta_t = ros::Time::now() - _time_start;

  all_legs_info.header.stamp = _zero_time + delta_t;

  for (size_t leg_num = 0; leg_num < 4; leg_num++)
  {
    all_legs_info.leg.at(leg_num).p_error.x =
      all_legs_info.leg.at(leg_num).p_des.x - all_legs_info.leg.at(leg_num).p_act.x;
    all_legs_info.leg.at(leg_num).p_error.y =
      all_legs_info.leg.at(leg_num).p_des.y - all_legs_info.leg.at(leg_num).p_act.y;
    all_legs_info.leg.at(leg_num).p_error.z =
      all_legs_info.leg.at(leg_num).p_des.z - all_legs_info.leg.at(leg_num).p_act.z;

    all_legs_info.leg.at(leg_num).v_error.x =
      all_legs_info.leg.at(leg_num).v_des.x - all_legs_info.leg.at(leg_num).v_act.x;
    all_legs_info.leg.at(leg_num).v_error.y =
      all_legs_info.leg.at(leg_num).v_des.y - all_legs_info.leg.at(leg_num).v_act.y;
    all_legs_info.leg.at(leg_num).v_error.z =
      all_legs_info.leg.at(leg_num).v_des.z - all_legs_info.leg.at(leg_num).v_act.z;
  }

  body_info.state_error.p.x = body_info.pos_des.x - body_info.pos_act.x;
  body_info.state_error.p.y = body_info.pos_des.y - body_info.pos_act.y;
  body_info.state_error.p.z = body_info.pos_des.z - body_info.pos_act.z;

  body_info.state_error.euler.x = body_info.euler_des.x - body_info.euler_act.x;
  body_info.state_error.euler.y = body_info.euler_des.y - body_info.euler_act.y;
  body_info.state_error.euler.z = body_info.euler_des.z - body_info.euler_act.z;

  body_info.state_error.v.linear.x = body_info.vel_des.linear.x - body_info.vel_act.linear.x;
  body_info.state_error.v.linear.y = body_info.vel_des.linear.y - body_info.vel_act.linear.y;
  body_info.state_error.v.linear.z = body_info.vel_des.linear.z - body_info.vel_act.linear.z;

  body_info.state_error.v.angular.x = body_info.vel_des.angular.x - body_info.vel_act.angular.x;
  body_info.state_error.v.angular.y = body_info.vel_des.angular.y - body_info.vel_act.angular.y;
  body_info.state_error.v.angular.z = body_info.vel_des.angular.z - body_info.vel_act.angular.z;

  _pub_all_legs_info.publish(all_legs_info);
  _pub_body_info.publish(body_info);
}

void Debug::updateVisualization()
{
  sensor_msgs::JointState msg;

  msg.header.stamp = ros::Time::now();

  msg.name.push_back("FR_hip_joint");
  msg.name.push_back("FR_thigh_joint");
  msg.name.push_back("FR_calf_joint");

  msg.name.push_back("FL_hip_joint");
  msg.name.push_back("FL_thigh_joint");
  msg.name.push_back("FL_calf_joint");

  msg.name.push_back("RR_hip_joint");
  msg.name.push_back("RR_thigh_joint");
  msg.name.push_back("RR_calf_joint");

  msg.name.push_back("RL_hip_joint");
  msg.name.push_back("RL_thigh_joint");
  msg.name.push_back("RL_calf_joint");

  for (size_t leg_num = 0; leg_num < 4; leg_num++)
  {
    for (uint8_t joint_num = 0; joint_num < 3; joint_num++)
    {
      msg.position.push_back(all_legs_info.leg.at(leg_num).joint.at(joint_num).q);
      msg.velocity.push_back(all_legs_info.leg.at(leg_num).joint.at(joint_num).dq);
    }
  }

  _pub_joint_states.publish(msg);
}

void Debug::tfPublish()
{
  geometry_msgs::TransformStamped odom_trans;

  odom_trans.header.stamp = ros::Time::now();
  odom_trans.header.frame_id = "odom";
  odom_trans.child_frame_id = "base";

  odom_trans.transform.translation.x = body_info.pos_act.x;
  odom_trans.transform.translation.y = body_info.pos_act.y;
  odom_trans.transform.translation.z = body_info.pos_act.z;

  geometry_msgs::Quaternion odom_quat;
  // TODO почему результаты естиматора приходится менять местами?
  odom_quat.x = body_info.quat_act.y;
  odom_quat.y = body_info.quat_act.z;
  odom_quat.z = body_info.quat_act.w;
  odom_quat.w = body_info.quat_act.x;
  odom_trans.transform.rotation = odom_quat;

  odom_broadcaster.sendTransform(odom_trans);

  geometry_msgs::TransformStamped odom_trans_world;

  odom_trans_world.header.stamp = odom_trans.header.stamp;
  odom_trans_world.header.frame_id = "world";
  odom_trans_world.child_frame_id = "odom";

  //  z_offset = _ground_trurh_pose.pose.pose.position.z - body_info.pos_act.z;
  odom_trans_world.transform.translation.z = z_offset;
  odom_trans_world.transform.rotation.w = 1.;

  world_odom_broadcaster.sendTransform(odom_trans_world);
}
