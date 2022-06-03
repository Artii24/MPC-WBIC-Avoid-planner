/*! @file PositionVelocityEstimator.cpp
 *  @brief All State Estimation Algorithms
 *
 *  This file will contain all state estimation algorithms.
 *  PositionVelocityEstimators should compute:
 *  - body position/velocity in world/body frames
 *  - foot positions/velocities in body/world frame
 */

#include "Controllers/MyPositionVelocityEstimator.h"

/*!
 * Initialize the state estimator
 */
template <typename T>
void KFPositionVelocityEstimator<T>::setup()
{
  T dt = 0.002;
  // T dt = this->_stateEstimatorData.parameters->controller_dt;
  _xhat.setZero();
  _ps.setZero();
  _vs.setZero();
  _A.setZero();
  _A.block(0, 0, 3, 3) = Eigen::Matrix<T, 3, 3>::Identity();
  _A.block(0, 3, 3, 3) = dt * Eigen::Matrix<T, 3, 3>::Identity();
  _A.block(3, 3, 3, 3) = Eigen::Matrix<T, 3, 3>::Identity();
  _A.block(6, 6, 12, 12) = Eigen::Matrix<T, 12, 12>::Identity();
  _B.setZero();
  _B.block(0, 0, 3, 3) = dt * dt * Eigen::Matrix<T, 3, 3>::Identity()/2;
  _B.block(3, 0, 3, 3) = dt * Eigen::Matrix<T, 3, 3>::Identity();
  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> C1(3, 6);
  C1 << Eigen::Matrix<T, 3, 3>::Identity(), Eigen::Matrix<T, 3, 3>::Zero();
  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> C2(3, 6);
  C2 << Eigen::Matrix<T, 3, 3>::Zero(), Eigen::Matrix<T, 3, 3>::Identity();
  _C.setZero();
  _C.block(0, 0, 3, 6) = C1;
  _C.block(3, 0, 3, 6) = C1;
  _C.block(6, 0, 3, 6) = C1;
  _C.block(9, 0, 3, 6) = C1;
  _C.block(0, 6, 12, 12) = T(-1) * Eigen::Matrix<T, 12, 12>::Identity();
  _C.block(12, 0, 3, 6) = C2;
  _C.block(15, 0, 3, 6) = C2;
  _C.block(18, 0, 3, 6) = C2;
  _C.block(21, 0, 3, 6) = C2;
  _C(27, 17) = T(1);
  _C(26, 14) = T(1);
  _C(25, 11) = T(1);
  _C(24, 8) = T(1);
  _P.setIdentity();
  _P = T(100) * _P;
  _Q0.setIdentity();
  _Q0.block(0, 0, 3, 3) = (dt / 20.f) * Eigen::Matrix<T, 3, 3>::Identity();
  _Q0.block(3, 3, 3, 3) =
      (dt * 9.8f / 20.f) * Eigen::Matrix<T, 3, 3>::Identity();
  _Q0.block(6, 6, 12, 12) = dt * Eigen::Matrix<T, 12, 12>::Identity();
  _R0.setIdentity();
  a_old = this->_stateEstimatorData.result->aWorld;
  a_filt <<0,0,0;
  da_filt <<0,0,0;
  a_world <<0,0,0;
  da_filt_prev << 0,0,0;
}

template <typename T>
KFPositionVelocityEstimator<T>::KFPositionVelocityEstimator() {}

/*!
 * Run state estimator
 */
template <typename T>
void KFPositionVelocityEstimator<T>::run()
{
  static uint16_t counter = 0;
    if (counter <= 2000)
  {
    counter++;
    a_old = this->_stateEstimatorData.result->aWorld;
  }
  else
  {
    T process_noise_pimu = this->_stateEstimatorData.parameters->imu_process_noise_position;
    T process_noise_vimu = this->_stateEstimatorData.parameters->imu_process_noise_velocity;
    T process_noise_pfoot = this->_stateEstimatorData.parameters->foot_process_noise_position;
    T sensor_noise_pimu_rel_foot = this->_stateEstimatorData.parameters->foot_sensor_noise_position;
    T sensor_noise_vimu_rel_foot = this->_stateEstimatorData.parameters->foot_sensor_noise_velocity;
    T sensor_noise_zfoot = this->_stateEstimatorData.parameters->foot_height_sensor_noise;
    // T contacts = this _>this->_stateEstimatorData->getContactSensorData()(leg_num);

    Eigen::Matrix<T, 18, 18> Q = Eigen::Matrix<T, 18, 18>::Identity();
    Q.block(0, 0, 3, 3) = _Q0.block(0, 0, 3, 3) * process_noise_pimu;
    Q.block(3, 3, 3, 3) = _Q0.block(3, 3, 3, 3) * process_noise_vimu;
    Q.block(6, 6, 12, 12) = _Q0.block(6, 6, 12, 12) * process_noise_pfoot;

    Eigen::Matrix<T, 28, 28> R = Eigen::Matrix<T, 28, 28>::Identity();
    R.block(0, 0, 12, 12) = _R0.block(0, 0, 12, 12) * sensor_noise_pimu_rel_foot;
    R.block(12, 12, 12, 12) = _R0.block(12, 12, 12, 12) * sensor_noise_vimu_rel_foot;
    R.block(24, 24, 4, 4) = _R0.block(24, 24, 4, 4) * sensor_noise_zfoot;

    int qindex = 0;
    int rindex1 = 0;
    int rindex2 = 0;
    int rindex3 = 0;

    Vec3<T> g(0, 0, T(-9.81));
    Mat3<T> Rbod = this->_stateEstimatorData.result->rBody.transpose();
    // in old code, Rbod * se_acc + g
    
    Vec3<T> a = this->_stateEstimatorData.result->aWorld;
    Vec3<T> w1(T(0.0001),T(0.0001),T(0.0001));
    Vec3<T> w2(T(0.00001),T(0.00001),T(0.00001));
    a_world += Vec3<T> (w1[0]*(a[0]-a_old[0])/(w1[0]+(a[0]-a_old[0])*(a[0]-a_old[0])),//(a[0]-a_old[0])*(a[0]-a_old[0])
                      w1[1]*(a[1]-a_old[1])/(w1[1]+(a[1]-a_old[1])*(a[1]-a_old[1])),//(a[1]-a_old[1])*(a[1]-a_old[1])
                      w1[2]*(a[2]-a_old[2])/(w1[2]+(a[2]-a_old[2])*(a[2]-a_old[2])));//(a[2]-a_old[2])*(a[2]-a_old[2])
    da_filt <<  T(1.0)*(a_world[0]-a_filt[0])*(a_world[0]-a_filt[0])*(a_world[0]-a_filt[0])*(a_world[0]-a_filt[0])*(a_world[0]-a_filt[0])/(w2[0] +(a_world[0]-a_filt[0])*(a_world[0]-a_filt[0])*(a_world[0]-a_filt[0])*(a_world[0]-a_filt[0])),
                T(1.0)*(a_world[1]-a_filt[1])*(a_world[1]-a_filt[1])*(a_world[1]-a_filt[1])*(a_world[1]-a_filt[1])*(a_world[1]-a_filt[1])/(w2[1] +(a_world[1]-a_filt[1])*(a_world[1]-a_filt[1])*(a_world[1]-a_filt[1])*(a_world[1]-a_filt[1])),
                T(1.0)*(a_world[2]-a_filt[2])*(a_world[2]-a_filt[2])*(a_world[2]-a_filt[2])*(a_world[2]-a_filt[2])*(a_world[2]-a_filt[2])/(w2[2] +(a_world[2]-a_filt[2])*(a_world[2]-a_filt[2])*(a_world[2]-a_filt[2])*(a_world[2]-a_filt[2]));
    a_old = a;
    a_filt = a_filt + 0.0*da_filt_prev + 1*da_filt;
    this->_stateEstimatorData.debug->body_info.test.x = a_filt[0];
    this->_stateEstimatorData.debug->body_info.test.y = a_filt[1];
    this->_stateEstimatorData.debug->body_info.test.z = a_filt[2];
    this->_stateEstimatorData.debug->body_info.test1.x = a_world[0];
    this->_stateEstimatorData.debug->body_info.test1.y = a_world[1];
    this->_stateEstimatorData.debug->body_info.test1.z = a_world[2];
    da_filt_prev = da_filt;
    
    // std::cout << "Error WORLD\n" << a[0] - a_filt[0] <<" " << a[1] - a_filt[1] <<" " << a[2] - a_filt[2] <<"\n";
    // std::cout << "A WORLD\n" << a[0] <<" " << a[1]  <<" " << a[2] <<"\n";
    // std::cout << "A filt WORLD\n" << a_filt[0] <<" " << a_filt[1] <<" " << a_filt[2] <<"\n";
    Vec4<T> pzs = Vec4<T>::Zero();
    Vec4<T> trusts = Vec4<T>::Zero();
    Vec3<T> p0, v0;
    Vec4<T> pzs0;
    p0 << _xhat[0], _xhat[1], _xhat[2];
    v0 << _xhat[3], _xhat[4], _xhat[5];
    pzs0 <<  _xhat[8], _xhat[11], _xhat[14], _xhat[17];

    for (int i = 0; i < 4; i++)
    {
      int i1 = 3 * i;
      Quadruped<T>& quadruped = *(this->_stateEstimatorData.legControllerData->quadruped);
      Vec3<T> ph = quadruped.getHipLocation(i); // hip positions relative to CoM
      // hw_i->leg_controller->leg_datas[i].p;
      Vec3<T> p_rel = ph + this->_stateEstimatorData.legControllerData[i].p; // Local frame distance from COM to leg(i)
      // hw_i->leg_controller->leg_datas[i].v;

      // Local frame velocity of leg(i) relative to COM
      Vec3<T> dp_rel = this->_stateEstimatorData.legControllerData[i].v;

      //Distance to leg(i) from Aligned with World frame body frame
      Vec3<T> p_f = Rbod * p_rel; 

      // World leg(i) velocity in Aligned with World frame body frame
      Vec3<T> dp_f = Rbod * (this->_stateEstimatorData.result->omegaBody.cross(p_rel) + dp_rel);
      

      qindex = 6 + i1;
      rindex1 = i1;
      rindex2 = 12 + i1;
      rindex3 = 24 + i;
      T k = 0;
      uint8_t a = (*this->_stateEstimatorData.contactSensor)[0];
      if (a) 
      {
        k=1;
      }
      T trust = T(1)*k;
      T phase = fmin(this->_stateEstimatorData.result->contactEstimate(i), T(1));
      //T trust_window = T(0.25);
      T trust_window = T(0.2);

      // trust from 0 to 1 when phase< trust_window and  phase > 1 - trust_window
      // In other phases trust = 1
      // That could happen because of the big noise/deviations during the contact moment
      if (phase < trust_window)
      {
        trust = phase / trust_window*k; // trust from 0 to 1
      }
      else if (phase > (T(1) - trust_window)) //  1 > (1- p)tr_w
      {
        trust = (T(1) - phase) / trust_window*k; // trust from 1 to 0
      }
      //T high_suspect_number(1000);
      T high_suspect_number(1000);

      // printf("Trust %d: %.3f\n", i, trust);
      Q.block(qindex, qindex, 3, 3) = (T(1) + (T(1) - trust) * high_suspect_number) * Q.block(qindex, qindex, 3, 3);
      R.block(rindex1, rindex1, 3, 3) = 1 * R.block(rindex1, rindex1, 3, 3);
      R.block(rindex2, rindex2, 3, 3) = (T(1) + (T(1) - trust) * high_suspect_number) * R.block(rindex2, rindex2, 3, 3);
      R(rindex3, rindex3) = (T(1) + (T(1) - trust) * high_suspect_number) * R(rindex3, rindex3);

      trusts(i) = trust;
    
      _ps.segment(i1, 3) = -p_f;
      _vs.segment(i1, 3) =  (1.0f - trust) *v0 + trust * (-dp_f);//
      pzs(i) = trust*pzs0(i)+ (1.0f - trust) *(p0(2) + p_f(2)) ;//
      std::cout <<"pzs("<< i << ") " <<pzs(i)<<" "<<"p0z "<< p0(2)<<" " << "phase " << phase << std::endl;
      std::cout <<"pzs("<< i << ") " <<pzs(i)<<" "<<"p0z "<< p0(2)<<" " << "phase " << phase << std::endl;
      std::cout <<"pzs("<< i << ") " <<pzs(i)<<" "<<"p0z "<< p0(2)<<" " << "phase " << phase << std::endl;

    }
    //std::cout <<"Trusts"<< " " << trusts(0)<< " "<< trusts(1)<< " "<< trusts(2)<< " "<< trusts(3)<< " "<< std::endl;
    
    Eigen::Matrix<T, 28, 1> y;
    y << _ps, _vs, pzs;
    _xhat = _A * _xhat + _B * a_filt;
    Eigen::Matrix<T, 18, 18> At = _A.transpose();
    Eigen::Matrix<T, 18, 18> Pm = _A * _P * At + Q;
    Eigen::Matrix<T, 18, 28> Ct = _C.transpose();
    Eigen::Matrix<T, 28, 1> yModel = _C * _xhat;
    Eigen::Matrix<T, 28, 1> ey = y - yModel;
    // std::cout << yModel[0] << " " << yModel[1] << " " << yModel[2] << std::endl; 
    Eigen::Matrix<T, 28, 28> S = _C * Pm * Ct + R;

    // todo compute LU only once
    // LU decompostion is a method to solve linear system of algebra equations, find inverse matrix 
    // and determinant
    // Ax=b or A.lu().solve(b) make lower-upper triangular matrix decomposition for A = LU
    // and the make a replacement z = Ux such that we have Lz=b
    // Example: 
    // A = [1 1;2 4] b = [5; 6]
    // L = [1 0;2 1] U = [1 1;0 2]
    // Z1 = 5 | 2*Z1+Z2 = 6 => Z2= -4
    // 2*X2 = -4 | X1 + X2 = 5 => X1 = 7
    // Check: 2*X1+4*X2 = 6 | 14-8=6 

    // A has to be invertiable.
    Eigen::Matrix<T, 28, 1> S_ey = S.lu().solve(ey);
    _xhat += Pm * Ct * S_ey;
    // std::cout << _xhat[2] << " ______________ "<<_xhat[2] << std::endl;
    // std::cout << _xhat[2] << " ______________ "<<_xhat[2] << std::endl;
    // std::cout << _xhat[2] << " ______________ "<<_xhat[2] << std::endl;
    Eigen::Matrix<T, 28, 18> S_C = S.lu().solve(_C);
    _P = (Eigen::Matrix<T, 18, 18>::Identity() - Pm * Ct * S_C) * Pm;

    Eigen::Matrix<T, 18, 18> Pt = _P.transpose();
    _P = (_P + Pt) / T(2);

    if (_P.block(0, 0, 2, 2).determinant() > T(0.000001))
    {
      _P.block(0, 2, 2, 16).setZero();
      _P.block(2, 0, 16, 2).setZero();
      _P.block(0, 0, 2, 2) /= T(10);
    }
    // std::cout << _xhat[2] << std::endl;
    //this->_stateEstimatorData.result->position = _xhat.block(0, 0, 3, 1);
    //this->_stateEstimatorData.result->vWorld = _xhat.block(3, 0, 3, 1);
    //this->_stateEstimatorData.result->vBody = this->_stateEstimatorData.result->rBody * this->_stateEstimatorData.result->vWorld;
  }
}
template class KFPositionVelocityEstimator<float>;
template class KFPositionVelocityEstimator<double>;

// /*!
//  * Run cheater estimator to copy cheater state into state estimate
//  */
// template <typename T>
// void CheaterPositionVelocityEstimator<T>::run()
// {
//   this->_stateEstimatorData.result->position = this->_stateEstimatorData.cheaterState->position.template cast<T>();
//   this->_stateEstimatorData.result->vWorld = this->_stateEstimatorData.result->rBody.transpose().template cast<T>() * this->_stateEstimatorData.cheaterState->vBody.template cast<T>();
//   this->_stateEstimatorData.result->vBody = this->_stateEstimatorData.cheaterState->vBody.template cast<T>();
// }

// template class CheaterPositionVelocityEstimator<float>;
// template class CheaterPositionVelocityEstimator<double>;
