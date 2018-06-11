#pragma once

#include <srrg_matchable/scene.h>

namespace srrg_g2o_simulator{

  struct CellPair{
    CellPair(const Eigen::Vector2i& first_cell_,
             const Eigen::Vector2i& second_cell_):
      first_cell(first_cell_),
      second_cell(second_cell_){}

    inline bool operator <(const CellPair &c) const {
      for(int i=0; i<2; ++i){
        if(first_cell[i] < c.first_cell[i])
          return true;
        if(first_cell[i] > c.first_cell[i])
          return false;
      }
      for(int i=0; i<2; ++i){
        if(second_cell[i] < c.second_cell[i])
          return true;
        if(second_cell[i] > c.second_cell[i])
          return false;
      }
      return false;
    }

    inline bool operator ==(const CellPair &c) const {
      for(int i=0; i<2; ++i)
        if(first_cell[i] != c.first_cell[i])
          return false;

      for(int i=0; i<2; ++i)
        if(second_cell[i] != c.second_cell[i])
          return false;

      return true;
    }

    Eigen::Vector2i first_cell;
    Eigen::Vector2i second_cell;
  };
  typedef std::map<CellPair,srrg_matchable::MatchablePtr> CellPairPlaneMap;
  typedef std::set<srrg_matchable::MatchablePtr> MatchablePtrSet;

  class Simulator{
  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

    inline void setResolution(float resolution_){_resolution=resolution_;}
    inline void setDimension(int d_){_dimension.x()=d_;_dimension.y()=d_;}
    inline void setNumPoses(int n_){_num_poses=n_;}
    inline void setNumPlanes(const int n_) {_num_planes = n_;}
    inline void setNumLines(const int n_) {_num_lines = n_;}
    inline void setNumPoints(const int n_) {_num_points = n_;}

    inline const srrg_matchable::MatchablePtrSet &matchables() const {return _set;}
    inline const std::vector<Eigen::Vector3f,
      Eigen::aligned_allocator<Eigen::Vector3f> >& robotTrajectory() const {return _trajectory;}

    void populateSet();

    void sbragation();

  protected:
    float _resolution;
    Eigen::Vector2i _dimension;
    int _num_poses;
    int _num_planes;
    int _num_lines;
    int _num_points;
    std::vector<Eigen::Vector3f, Eigen::aligned_allocator<Eigen::Vector3f> > _trajectory;

    srrg_matchable::MatchablePtrSet _set;
    CellPairPlaneMap _map;

  private:
    void computeRotationMatrix(Eigen::Matrix3f& rotation_matrix,
                               const Eigen::Vector3f& direction);
  };
}
