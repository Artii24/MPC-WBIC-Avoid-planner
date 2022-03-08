/*!
 * @file RobotController.h
 * @brief Parent class of user robot controllers.
 * This is an interface between the control code and the common hardware code
 */

#ifndef ROBOT_CONTROLLER_H
#define ROBOT_CONTROLLER_H

#include "Controllers/LegController.h"
#include "Controllers/StateEstimatorContainer.h"
#include "Dynamics/FloatingBaseModel.h"
// #include "SimUtilities/VisualizationData.h"

/*!
 * Parent class of user robot controllers
 */
class RobotController
{

public:
  RobotController() {}
  virtual ~RobotController() {}

  virtual void initializeController() = 0;
  /**
 * Called one time every control loop 
 */
  virtual void runController() = 0;
  virtual ControlParameters* getUserControlParameters() = 0;
  virtual void Estop() {}

  Quadruped<float>* _quadruped = nullptr;
  FloatingBaseModel<float>* _model = nullptr;
  LegController<float>* _legController = nullptr;
  StateEstimatorContainer<float>* _stateEstimator = nullptr;
  StateEstimate<float>* _stateEstimate = nullptr;
  GamepadCommand* _driverCommand = nullptr;
  RobotControlParameters* _controlParameters = nullptr;
  DesiredStateCommand<float>* _desiredStateCommand = nullptr;

  // VisualizationData* _visualizationData = nullptr;
  RobotType _robotType;
};

#endif
