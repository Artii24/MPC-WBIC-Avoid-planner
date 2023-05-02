#include "WBC_Ctrl.hpp"
#include "FloatingBaseModel.h"
#include "WBC/WBC.hpp"
#include <Utilities/Timer.h>
#include <Utilities/Utilities_print.h>

template<typename T>
WBC_Ctrl<T>::WBC_Ctrl(FloatingBaseModel<T> model)
  : _full_config(cheetah::num_act_joint + 7),
    _tau_ff(cheetah::num_act_joint),
    _des_jpos(cheetah::num_act_joint),
    _des_jvel(cheetah::num_act_joint)
{
  _iter = 0;
  _full_config.setZero();

  _model = model;
  _kin_wbc = new KinWBC<T>(cheetah::dim_config);

  _wbic = new WBIC<T>(cheetah::dim_config, &(_contact_list), &(_task_list));
  _wbic_data = new WBIC_ExtraData<T>();

  _wbic_data->_W_floating = DVec<T>::Constant(6, 0.1);

  _wbic_data->_W_rf = DVec<T>::Constant(12, 1.);

  _Kp_joint = DVec<T>::Constant(cheetah::num_leg_joint, 5.);
  _Kd_joint = DVec<T>::Constant(cheetah::num_leg_joint, 1.5);

  _state.q = DVec<T>::Zero(cheetah::num_act_joint);
  _state.qd = DVec<T>::Zero(cheetah::num_act_joint);
}

template<typename T>
WBC_Ctrl<T>::~WBC_Ctrl()
{
  delete _kin_wbc;
  delete _wbic;
  delete _wbic_data;

  typename std::vector<Task<T>*>::iterator iter = _task_list.begin();
  while (iter < _task_list.end())
  {
    delete (*iter);
    ++iter;
  }
  _task_list.clear();

  typename std::vector<ContactSpec<T>*>::iterator iter2 = _contact_list.begin();
  while (iter2 < _contact_list.end())
  {
    delete (*iter2);
    ++iter2;
  }
  _contact_list.clear();
}

template<typename T>
void WBC_Ctrl<T>::_ComputeWBC()
{
  // TEST
  _kin_wbc->FindConfiguration(_full_config, _task_list, _contact_list, _des_jpos, _des_jvel);

  // WBIC
  _wbic->UpdateSetting(_A, _Ainv, _coriolis, _grav);
  _wbic->MakeTorque(_tau_ff, _wbic_data);
}

template<typename T>
void WBC_Ctrl<T>::run(void* input, ControlFSMData<T>& data)
{
  ++_iter;

  // Update Model
  _UpdateModel(data.stateEstimator->getResult(), data.legController->datas);

  // Task & Contact Update
  _ContactTaskUpdate(input, data);

  // WBC Computation
  _ComputeWBC();

  // Update Leg Command
  Vec4<T> contact = data._stateEstimator->getResult().contactEstimate;
  int shift = 0;
  for (size_t leg(0); leg < cheetah::num_leg; ++leg)
  {
    // std::cout<<contact(leg)<<endl;
    if (shift<2)
    {

    
    if (contact(leg) > 0.01)
    {

      data.debug->all_legs_info.leg[leg].wbc_force.x = -_wbic_data->_Fr(cheetah::num_leg_joint * shift );
      data.debug->all_legs_info.leg[leg].wbc_force.y = -_wbic_data->_Fr(cheetah::num_leg_joint * shift +1);
      data.debug->all_legs_info.leg[leg].wbc_force.z = -_wbic_data->_Fr(cheetah::num_leg_joint * shift +2);
      shift = shift +1;
      // std::cout<<shift<<endl;
    }
    else
    {
      data.debug->all_legs_info.leg[leg].wbc_force.x = 0;
      data.debug->all_legs_info.leg[leg].wbc_force.y = 0;
      data.debug->all_legs_info.leg[leg].wbc_force.z = 0;
    }
    }
    // data.debug->all_legs_info.leg[i].wbc_force.x = _wbic_data->_opt_result(3*i);
    // data.debug->all_legs_info.leg[i].wbc_force.y = _wbic_data->_opt_result(3*i+1);
    // data.debug->all_legs_info.leg[i].wbc_force.z = _wbic_data->_opt_result(3*i+2);
    // data.debug->all_legs_info.leg[i].wbc_force.w = std::sqrt(_wbic_data->_opt_result.block(3*i,0,3,1).dot(_wbic_data->_opt_result.block(3*i,0,3,1)));
    // data.debug->all_legs_info.leg[i].wbc_force.x = _wbic_data->_Fr(3);
  }
  
  _UpdateLegCMD(data);
}

template<typename T>
void WBC_Ctrl<T>::_UpdateLegCMD(ControlFSMData<T>& data)
{
  LegControllerCommand<T>* cmd = data._legController->commands;
  Vec4<T> contact = data._stateEstimator->getResult().contactEstimate;
  int shift = 0;
  for (size_t leg(0); leg < cheetah::num_leg; ++leg)
  {
    cmd[leg].zero();
    // if (contact(leg) == 1)
    // {
    //   data.debug->all_legs_info.leg[leg].wbc_force.x = -_wbic_data->_Fr(cheetah::num_leg_joint * shift );
    //   data.debug->all_legs_info.leg[leg].wbc_force.y = -_wbic_data->_Fr(cheetah::num_leg_joint * shift +1);
    //   data.debug->all_legs_info.leg[leg].wbc_force.z = -_wbic_data->_Fr(cheetah::num_leg_joint * shift +2);
    //   shift = shift +1;
    // }
    for (size_t jidx(0); jidx < cheetah::num_leg_joint; ++jidx)
    {
      cmd[leg].tauFeedForward[jidx] = _tau_ff[cheetah::num_leg_joint * leg + jidx];


      cmd[leg].qDes[jidx] = _des_jpos[cheetah::num_leg_joint * leg + jidx];
      cmd[leg].qdDes[jidx] = _des_jvel[cheetah::num_leg_joint * leg + jidx];

      cmd[leg].kdJoint(jidx, jidx) = _Kd_joint[jidx];
      cmd[leg].kpJoint(jidx, jidx) = _Kp_joint[jidx];
    }
  }

  // Knee joint non flip barrier
  for (size_t leg(0); leg < 4; ++leg)
  {
    if (cmd[leg].qDes[2] < 0.3)
    {
      cmd[leg].qDes[2] = 0.3;
    }
    if (data.legController->datas[leg].q[2] < 0.3)
    {
      T knee_pos = data.legController->datas[leg].q[2];
      cmd[leg].tauFeedForward[2] = 1. / (knee_pos * knee_pos + 0.02);
    }
  }
}

template<typename T>
void WBC_Ctrl<T>::_UpdateModel(const StateEstimate<T>& state_est, const LegControllerData<T>* leg_data)
{
  _state.bodyOrientation = state_est.orientation;
  _state.bodyPosition = state_est.position;

  for (size_t i(0); i < 3; ++i)
  {
    _state.bodyVelocity[i] = state_est.omegaBody[i];
    _state.bodyVelocity[i + 3] = state_est.vBody[i];

    for (size_t leg(0); leg < 4; ++leg)
    {
      _state.q[3 * leg + i] = leg_data[leg].q[i];
      _state.qd[3 * leg + i] = leg_data[leg].qd[i];

      _full_config[3 * leg + i + 6] = _state.q[3 * leg + i];
    }
  }
  _model.setState(_state);

  _model.contactJacobians();
  _model.massMatrix();
  _model.generalizedGravityForce();
  _model.generalizedCoriolisForce();

  _A.setZero(18, 18);
  _A = _model.getMassMatrix();
  _grav = _model.getGravityForce();
  _coriolis = _model.getCoriolisForce();
  _Ainv = _A.inverse();
}

template class WBC_Ctrl<float>;
template class WBC_Ctrl<double>;
