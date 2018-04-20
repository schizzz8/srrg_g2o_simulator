add_library(srrg_g2o_simulator_library SHARED
  simulator.cpp simulator.h
)

target_link_libraries(srrg_g2o_simulator_library
  ${catkin_LIBRARIES}
)
