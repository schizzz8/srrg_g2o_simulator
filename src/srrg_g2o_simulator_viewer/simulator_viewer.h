#pragma once

#include <QKeyEvent>
#include <QGLViewer/qglviewer.h>

#include <srrg_matchable/scene.h>


namespace srrg_g2o_simulator {

  class SimulatorViewer: public QGLViewer {
  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    //! ctor
    SimulatorViewer(QWidget* parent = 0);

    //! init method, opens the gl viewport and sets the key bindings
    void init();

    //! callback invoked by the application on new key event. It saves the last event in
    //! a member variable
    virtual void keyPressEvent(QKeyEvent *e);

    //! returns the last key pressed since invoking keyEventProcessed();
    QKeyEvent* lastKeyEvent();

    //! call this to clear the events, after processing them
    void keyEventProcessed();

    virtual void draw();

    inline void setScene(srrg_matchable::Scene* scene_){_scene=scene_;}
    inline void setRobotTrajectory(const std::vector<Eigen::Vector3f,
                                   Eigen::aligned_allocator<Eigen::Vector3f> >& trajectory_) {
      _robot_trajectory=trajectory_;
    }

  protected:

    QKeyEvent _last_key_event;
    bool _last_key_event_processed;

    srrg_matchable::Scene* _scene;
    std::vector<Eigen::Vector3f, Eigen::aligned_allocator<Eigen::Vector3f> > _robot_trajectory;

  };
}
