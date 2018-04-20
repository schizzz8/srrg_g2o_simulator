add_executable(srrg_g2o_simulator_app srrg_g2o_simulator_app.cpp)

target_link_libraries(srrg_g2o_simulator_app
  srrg_g2o_simulator_library
  srrg_g2o_simulator_viewer_library
  srrg_matchable_library
  ${OpenCV_LIBS}
  ${catkin_LIBRARIES}
)
