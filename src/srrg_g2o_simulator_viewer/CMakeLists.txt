add_library(srrg_g2o_simulator_viewer_library SHARED 
  simulator_viewer.cpp simulator_viewer.h
)

target_link_libraries(srrg_g2o_simulator_viewer_library
  ${SRRG_QT_LIBRARIES}
  ${QGLVIEWER_LIBRARY}
  ${catkin_LIBRARIES}
)
