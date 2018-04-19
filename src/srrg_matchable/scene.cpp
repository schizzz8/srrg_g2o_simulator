#include <iostream>
#include <GL/gl.h>
#include "scene.h"

namespace srrg_matchable {

  Scene::~Scene() {
    clear();
  }
  
  void Scene::draw(bool show_direction_matrix,
                     const Eigen::Vector3f& base_color) {
    int i=0;
    for(MatchablePtr m : *this) {
      if(!m)
        std::cerr << "[SRRG_MATCHABLE][SCENE][DRAW]: null pointer" << i << std::endl;
      m->draw(show_direction_matrix, base_color);
      i++;
    }
  }
  
  void Scene::drawWithNames(int name_offset,
                              bool show_direction_matrix,
                              const Eigen::Vector3f& base_color) {
    int name_idx = name_offset;
    for(MatchablePtr m : *this) {
      glPushName(name_idx);
      m->draw(show_direction_matrix, base_color);
      glPopName();
      ++name_idx;
    }    
  }  

  void Scene::transformInPlace(const Eigen::Isometry3f& iso) const {
    for(MatchablePtr m : *this)
      m->transformInPlace(iso);    
  }
  
  Scene Scene::transform(const Eigen::Isometry3f& iso) const {
    Scene s2 = *this;
    s2.transformInPlace(iso);
    return s2;    
  }
  
  Scene& Scene::operator=(const Scene& s2) {
    clear();
    if(! s2.size())
      return *this;
    resize(s2.size());
    for(size_t i = 0; i < size(); ++i)
      at(i) = s2.at(i); // no clone since the pointer is shared
    return *this;    
  }


  Scene* Scene::clone() const {
    Scene* other_scene = new Scene();
    other_scene->resize(size());
    for(size_t i = 0; i < size(); ++i)
      other_scene->at(i) = at(i); // here we don't clone since the pointer is shared
    return other_scene;
  }

  void Scene::add(const Scene* other){
    if (! other)
      return;
    size_t old_size=size();
    resize(old_size+other->size());
    size_t k=old_size;
    for(const MatchablePtr m: *other){
      at(k)=m; // here we don't clone since the pointer is shared
      k++;
    }
  }
  
}
