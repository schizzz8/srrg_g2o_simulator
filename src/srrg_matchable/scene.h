#pragma once
#include "matchable.h"
#include <vector>

namespace srrg_matchable{
  class Scene : public std::vector<MatchablePtr> {
  public:
    virtual ~Scene();
    Scene* clone() const;
    Scene& operator=(const Scene& s2);
    void draw(bool show_direction_matrix=false,
              const Eigen::Vector3f& base_color = Eigen::Vector3f::Zero());
    void drawWithNames(int name_offset = 0, bool show_direction_matrix = false,
                       const Eigen::Vector3f& base_color = Eigen::Vector3f::Zero());
    Scene transform(const Eigen::Isometry3f& iso) const;
    void transformInPlace(const Eigen::Isometry3f& iso) const;
    void add(const Scene* other);
    
  };
}
