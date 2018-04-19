add_library(srrg_matchable_library SHARED
  matchable.cpp matchable.h
  scene.cpp scene.h
  )

target_link_libraries(srrg_matchable_library
  ${catkin_LIBRARIES}
)
