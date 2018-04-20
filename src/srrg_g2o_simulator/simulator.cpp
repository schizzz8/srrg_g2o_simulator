#include "simulator.h"

namespace srrg_g2o_simulator{

  using namespace srrg_matchable;

  void Simulator::populateSet(){

    //insert points
    for(size_t j=0; j<_dimension.y(); ++j)
      for(size_t i=0; i<_dimension.x(); ++i){
        const Eigen::Vector3f p(i*_resolution, j*_resolution, 0.0f);
        MatchablePtr point_matchable(new Matchable(Matchable::Point,p));
        _set.insert(point_matchable);
      }

    //insert lines along x axes
    for(size_t j=0; j<_dimension.y(); ++j) //
      for(size_t i=0; i<_dimension.x()-1; ++i){
        const Eigen::Vector3f p(i*_resolution, j*_resolution, 0.0f);
        const Eigen::Vector3f n(1.0f,0.0f,0.0f);
        MatchablePtr line_matchable(new Matchable(Matchable::Line,p,n,Eigen::Vector2f(_resolution,0.0f)));
        Eigen::Matrix3f R;
        computeRotationMatrix(R,n);
        line_matchable->setRotationMatrix(R);
        _set.insert(line_matchable);
      }

    //insert lines along y axes
    for(size_t i=0; i<_dimension.x(); ++i)
      for(size_t j=0; j<_dimension.y()-1; ++j){
        const Eigen::Vector3f p(i*_resolution, j*_resolution, 0.0f);
        const Eigen::Vector3f n(0.0f,1.0f,0.0f);
        MatchablePtr line_matchable(new Matchable(Matchable::Line,p,n,Eigen::Vector2f(_resolution,0.0f)));
        Eigen::Matrix3f R;
        computeRotationMatrix(R,n);
        line_matchable->setRotationMatrix(R);
        _set.insert(line_matchable);
      }

    //insert lines along z axes
    for(size_t j=0; j<_dimension.y(); ++j)
      for(size_t i=0; i<_dimension.x(); ++i){
        const Eigen::Vector3f p(i*_resolution, j*_resolution, 0.0f);
        const Eigen::Vector3f n(0.0f,0.0f,1.0f);
        MatchablePtr line_matchable(new Matchable(Matchable::Line,p,n,Eigen::Vector2f(_resolution,0.0f)));
        Eigen::Matrix3f R;
        computeRotationMatrix(R,n);
        line_matchable->setRotationMatrix(R);
        _set.insert(line_matchable);
      }

    //insert planes along x axes
    for(size_t j=0; j<_dimension.y()-1; ++j)
      for(size_t i=0; i<_dimension.x(); ++i){
        const Eigen::Vector2i cell(i,j);
        const Eigen::Vector3f p(cell.x()*_resolution, cell.y()*_resolution+_resolution/2.0f, _resolution/2.0f);
        const Eigen::Vector3f n(1.0f,0.0f,0.0f);
        MatchablePtr plane_matchable(new Matchable(Matchable::Plane,p,n,Eigen::Vector2f(_resolution/2.0f,_resolution/2.0f)));
        Eigen::Matrix3f R;
        computeRotationMatrix(R,n);
        plane_matchable->setRotationMatrix(R);
        _set.insert(plane_matchable);

        if(i >= 0 && i < _dimension.x()-1){
          const Eigen::Vector2i left_cell(i-1,j);
          CellPair left_pair(left_cell,cell);
          _map[left_pair] = plane_matchable;
          CellPair right_pair(cell,left_cell);
          _map[right_pair] = plane_matchable;
        }
      }

    //insert planes along y axes
    for(size_t i=0; i<_dimension.x()-1; ++i)
      for(size_t j=0; j<_dimension.y(); ++j){
        const Eigen::Vector2i cell(i,j);
        const Eigen::Vector3f p(cell.x()*_resolution+_resolution/2.0f, cell.y()*_resolution, _resolution/2.0f);
        const Eigen::Vector3f n(0.0f,1.0f,0.0f);
        MatchablePtr plane_matchable(new Matchable(Matchable::Plane,p,n,Eigen::Vector2f(_resolution/2.0f,_resolution/2.0f)));
        Eigen::Matrix3f R;
        computeRotationMatrix(R,n);
        plane_matchable->setRotationMatrix(R);
        _set.insert(plane_matchable);

        if(j >= 0 && j < _dimension.y()-1){
          const Eigen::Vector2i down_cell(i,j-1);
          CellPair down_pair(down_cell,cell);
          _map[down_pair] = plane_matchable;
          CellPair up_pair(cell,down_cell);
          _map[up_pair] = plane_matchable;
        }
      }

    //insert planes along z axes
    for(size_t j=0; j<_dimension.y()-1; ++j)
      for(size_t i=0; i<_dimension.x()-1; ++i){
        const Eigen::Vector3f p(i*_resolution+_resolution/2.0f, j*_resolution+_resolution/2.0f, 0.0f);
        const Eigen::Vector3f n(0.0f,0.0f,1.0f);
        MatchablePtr plane_matchable(new Matchable(Matchable::Plane,p,n,Eigen::Vector2f(_resolution/2.0f,_resolution/2.0f)));
        Eigen::Matrix3f R;
        computeRotationMatrix(R,n);
        plane_matchable->setRotationMatrix(R);
        _set.insert(plane_matchable);

      }
    std::cerr << "Scene has " << _set.size() << " matchables before sbraco" << std::endl;
  }

  void Simulator::sbragation(){
    std::cerr << "sbraco..." << std::endl;
    int count=0;
    bool continue_=true;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0, 1);

    Eigen::Vector3f position(_dimension.x()/2-1,_dimension.y()/2-1,0.0f);

    while(continue_){

      //sample new position
      Eigen::Vector3f increment = Eigen::Vector3f::Zero();
      float n = dis(gen);

      //go forward
      if(n < 0.5f){
        increment.x() = round(cos(position.z()));
        increment.y() = round(sin(position.z()));
      }
      //go left
      if(n >= 0.5f && n < 0.75f){
        increment.x() = round(-sin(position.z()));
        increment.y() = round(cos(position.z()));
        increment.z() = M_PI/2.0f;
      }
      //go right
      if(n >= 0.75f){
        increment.x() = round(sin(position.z()));
        increment.y() = round(-cos(position.z()));
        increment.z() = -M_PI/2.0f;
      }
      Eigen::Vector3f new_position = position+increment;

      //check if new position is out of grid
      if(new_position.x() < 0.0f || new_position.x() >= (float)(_dimension.x()-1) ||
         new_position.y() < 0.0f || new_position.y() >= (float)(_dimension.y()-1)){
        continue;
      }

      //sbraco
      const Eigen::Vector2i p = position.head(2).cast<int>();
      const Eigen::Vector2i np = new_position.head(2).cast<int>();
      CellPair pair(p,np);
      CellPairPlaneMap::iterator it = _map.find(pair);

      if(it != _map.end()){
        MatchablePtr m = it->second;
        MatchablePtrSet::iterator jt = _set.find(m);
        if(jt != _set.end())
          _set.erase(jt);
      } else {
        std::cerr << pair.first_cell.transpose() << " - " << pair.second_cell.transpose() << " not found!" << std::endl;
      }

      position = new_position;
      count++;

      if(count == _num_poses)
        continue_=false;

    }
  }

  void Simulator::computeRotationMatrix(Eigen::Matrix3f &rotation_matrix, const Eigen::Vector3f &direction){
    float d = sqrt(direction.x()*direction.x() + direction.y()*direction.y());

    const float& dirx = direction.x();
    const float& diry = direction.y();
    const float& dirz = direction.z();

    if(d > std::numeric_limits<float>::min()) {
      rotation_matrix << diry/d,     dirx*dirz/d,     dirx,
          -dirx/d,   diry*dirz/d, diry,
          0,          -d,        dirz;
    } else {
      rotation_matrix.setIdentity();
    }
  }

}
