#include "debug.hpp"

Debug::Debug(ros::Time time_start)
    : _zero_time(0), _time_start(time_start)
{
  z_offset = 0.f;
  _init();
}

void Debug::_init()
{
  _initPublishers();

  body_info.quat_act.w = 1.0;
}

void Debug::_initPublishers()
{
  _pub_joint_states = _nh.advertise<sensor_msgs::JointState>("/joint_states", 1);
  _pub_all_legs_info = _nh.advertise<unitree_legged_msgs::AllLegsInfo>("/all_legs_info", 1);
  _pub_body_info = _nh.advertise<unitree_legged_msgs::BodyInfo>("/body_info", 1);

  _pub_visual_last_p_stance = _nh.advertise<visualization_msgs::Marker>("/visual/last_p_stance", 1);
  _pub_estimated_stance_plane = _nh.advertise<visualization_msgs::Marker>("/visual/estimated_stance_plane", 1);

#ifdef PUB_IMU_AND_ODOM
  _pub_odom = _nh.advertise<nav_msgs::Odometry>("/odom", 1);
  _pub_imu = _nh.advertise<sensor_msgs::Imu>("/imu", 1);
#endif
}

void Debug::updatePlot()
{
  ros::Duration delta_t = ros::Time::now() - _time_start;

  // all_legs_info.header.stamp = _zero_time + delta_t;
  // body_info.header.stamp = _zero_time + delta_t;
  all_legs_info.header.stamp = ros::Time::now();
  body_info.header.stamp = ros::Time::now();

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

#ifdef PUB_IMU_AND_ODOM
  nav_msgs::Odometry odom;

  odom.header.stamp = ros::Time::now();
  odom.header.frame_id = "odom";
  odom.child_frame_id = "base";

  odom.pose.pose.position = body_info.pos_act;

  geometry_msgs::Quaternion odom_quat;
  odom_quat.x = body_info.quat_act.y;
  odom_quat.y = body_info.quat_act.z;
  odom_quat.z = body_info.quat_act.w;
  odom_quat.w = body_info.quat_act.x;
  odom.pose.pose.orientation = odom_quat;
  odom.twist.twist.linear = body_info.vel_act.linear;
  odom.twist.twist.angular = body_info.vel_act.angular;

  imu.header.frame_id = "imu_link";
  imu.header.stamp = ros::Time::now();

  imu.orientation = odom_quat;

  _pub_odom.publish(odom);
  _pub_imu.publish(imu);
#endif
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

  _drawLastStancePoints();
  _drawEstimatedStancePLane();
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

  z_offset = ground_truth_odom.pose.pose.position.z - body_info.pos_act.z;
  odom_trans_world.transform.translation.z = z_offset;
  odom_trans_world.transform.rotation.w = 1.;

  world_odom_broadcaster.sendTransform(odom_trans_world);
}

Vec3<float> Debug::_getHipLocation(uint8_t leg_num)
{
  Vec3<float> result(0, 0, 0);

  float x = 0.1805;
  float y = 0.047;

  switch (leg_num)
  {
    // FR
  case 0:
    result(0) = x;
    result(1) = -y;
    break;

    // FL
  case 1:
    result(0) = x;
    result(1) = y;
    break;

    // BR
  case 2:
    result(0) = -x;
    result(1) = -y;
    break;

    //BL
  case 3:
    result(0) = -x;
    result(1) = y;
    break;
  }

  return result;
}

void Debug::_drawLastStancePoints()
{
  visualization_msgs::Marker marker;

  std::string name = "odom";

  marker.header.frame_id = name;
  marker.header.stamp = ros::Time::now();
  marker.id = 0;
  marker.type = visualization_msgs::Marker::SPHERE_LIST;
  marker.action = visualization_msgs::Marker::ADD;
  // pose and orientation must be zero, except orientation.w = 1
  marker.pose.position.x = 0;
  marker.pose.position.y = 0;
  marker.pose.position.z = 0;
  marker.pose.orientation.x = 0.0;
  marker.pose.orientation.y = 0.0;
  marker.pose.orientation.z = 0.0;
  marker.pose.orientation.w = 1.0;
  marker.scale.x = 0.05;
  marker.scale.y = 0.05;
  marker.scale.z = 0.05;
  marker.color.a = 1.0;
  marker.color.r = 0.0;
  marker.color.g = 1.0;
  marker.color.b = 1.0;

  geometry_msgs::Point p0, p1, p2, p3;
  p0 = last_p_stance[0]; // FL point
  p1 = last_p_stance[1]; // FR point
  p2 = last_p_stance[2]; // BR point
  p3 = last_p_stance[3]; // BL point
  marker.points.clear();
  marker.points.push_back(p0);
  marker.points.push_back(p1);
  marker.points.push_back(p3);
  marker.points.push_back(p2);

  _pub_visual_last_p_stance.publish(marker);
}

void Debug::_drawEstimatedStancePLane()
{
  visualization_msgs::Marker marker;

  std::string name = "odom";

  tf::Quaternion quat;
  quat.setRPY(body_info.euler_act.x, body_info.euler_act.y, body_info.euler_act.z);

  marker.header.frame_id = name;
  marker.header.stamp = ros::Time::now();
  marker.id = 0;
  marker.type = visualization_msgs::Marker::CUBE;
  marker.action = visualization_msgs::Marker::ADD;
  // pose and orientation must be zero, except orientation.w = 1
  marker.pose.position = body_info.pos_act;
  marker.pose.position.z = 0;
  marker.pose.orientation.x = quat.x();
  marker.pose.orientation.y = quat.y();
  marker.pose.orientation.z = quat.z();
  marker.pose.orientation.w = quat.w();
  marker.scale.x = 0.5;
  marker.scale.y = 0.5;
  marker.scale.z = 0.00001;
  marker.color.a = 1.0;
  marker.color.r = 0.0;
  marker.color.g = 1.0;
  marker.color.b = 0.0;

  _pub_estimated_stance_plane.publish(marker);
}