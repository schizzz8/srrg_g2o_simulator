#include "simulator.h"

namespace srrg_g2o_simulator{

  using namespace srrg_matchable;

  void Simulator::populateSet(){
    //ia random number generator
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> pos_distribution(0, 1);
    std::uniform_real_distribution<> axis_distribution(0,2.5);

    //ia insert points
    for (size_t i = 0; i < _num_points; ++i) {
      float x = pos_distribution(gen)*_dimension.x();
      float y = pos_distribution(gen)*_dimension.y();
      float z = pos_distribution(gen);

      MatchablePtr point_matchable(new Matchable(Matchable::Point,
                                                 Eigen::Vector3f(x,y,z)));
      _set.insert(point_matchable);
    }


    //ia new lines generation
    //ia TODO: probabilistic axis selection
    //ia along x
    for (size_t i = 0; i < _num_lines/3; ++i) {
      float x = pos_distribution(gen)*_dimension.x();
      float y = pos_distribution(gen)*_dimension.y();
      float z = pos_distribution(gen);

      MatchablePtr line_matchable(new Matchable(Matchable::Line,
                                                Eigen::Vector3f(x,y,z),
                                                Eigen::Vector3f::UnitX(),
                                                Eigen::Vector2f(_resolution,0.0f)));
      Eigen::Matrix3f R;
      computeRotationMatrix(R,Eigen::Vector3f::UnitX());
      line_matchable->setRotationMatrix(R);
      _set.insert(line_matchable);
    }

    std::cerr << std::endl;
    
    //ia along y
    for (size_t i = 0; i < _num_lines/3; ++i) {
      float x = pos_distribution(gen)*_dimension.x();
      float y = pos_distribution(gen)*_dimension.y();
      float z = pos_distribution(gen);
        
      MatchablePtr line_matchable(new Matchable(Matchable::Line,
                                                Eigen::Vector3f(x,y,z),
                                                Eigen::Vector3f::UnitY(),
                                                Eigen::Vector2f(_resolution,0.0f)));
      Eigen::Matrix3f R;
      computeRotationMatrix(R,Eigen::Vector3f::UnitY());
      line_matchable->setRotationMatrix(R);
      _set.insert(line_matchable);
    }

    std::cerr << std::endl;

    //ia along z
    for (size_t i = 0; i < _num_lines/3; ++i) {
      float x = pos_distribution(gen)*_dimension.x();
      float y = pos_distribution(gen)*_dimension.y();
      float z = pos_distribution(gen);
      
      MatchablePtr line_matchable(new Matchable(Matchable::Line,
                                                Eigen::Vector3f(x,y,z),
                                                Eigen::Vector3f::UnitZ(),
                                                Eigen::Vector2f(_resolution,0.0f)));
      Eigen::Matrix3f R;
      computeRotationMatrix(R,Eigen::Vector3f::UnitZ());
      line_matchable->setRotationMatrix(R);
      _set.insert(line_matchable);
    }

    std::cerr << std::endl;

    
    //ia new planes generation
    int i = 0;
    while (i < _num_planes) {
      float x = 0;
      float y = 0;
      float z = 0;

      //ia grid coords
      int r = 0, c = 0;
      
      Eigen::Vector3f n=Eigen::Vector3f::Zero();
      Eigen::Vector2i prev_cell = Eigen::Vector2i::Zero();
      Eigen::Vector2i curr_cell = Eigen::Vector2i::Zero();

      const int axis_selector = round(axis_distribution(gen));
      switch (axis_selector) {
      case (0):
        {
          n = Eigen::Vector3f::UnitX();
          r = round(pos_distribution(gen)*_dimension.x());
          c = round(pos_distribution(gen)*_dimension.y());

          x = r * _resolution;
          y = c * _resolution + _resolution/2.0;
          
          z = pos_distribution(gen);

          prev_cell.x() = r - 1;
          prev_cell.y() = c;
          curr_cell.x() = r;
          curr_cell.y() = c;
          break;
        }
      case (1):
        {
          n = Eigen::Vector3f::UnitY();
          r = round(pos_distribution(gen)*_dimension.x());
          c = round(pos_distribution(gen)*_dimension.y());

          x = r * _resolution + _resolution/2.0;
          y = c * _resolution;
          
          z = pos_distribution(gen);

          prev_cell.x() = r;
          prev_cell.y() = c - 1;
          curr_cell.x() = r;
          curr_cell.y() = c;
          break;
        }
      case (2):
        {
          n = Eigen::Vector3f::UnitZ();
          x = round(pos_distribution(gen)*_dimension.x()) + _resolution/2.0f;
          y = round(pos_distribution(gen)*_dimension.y()) + _resolution/2.0f;
          break;
        }
      default:
        break;
      }

      if (x - _resolution/2 < 0 || x + _resolution/2 >= _dimension.x() ||
          y - _resolution/2 < 0 || y + _resolution/2 >= _dimension.y() ||
          z - _resolution/2 < 0 && axis_selector!= 2) {
        continue;
      }
      
      Eigen::Vector3f p(x,y,z);
      MatchablePtr plane_matchable(new Matchable(Matchable::Plane,
                                                 p,
                                                 n,
                                                 Eigen::Vector2f(_resolution/2.0f,
                                                                 _resolution/2.0f)));
      Eigen::Matrix3f R;
      computeRotationMatrix(R,n);
      plane_matchable->setRotationMatrix(R);
      _set.insert(plane_matchable);


      //ia wall
      if (axis_selector!=2) {
        CellPair pair0(prev_cell,curr_cell);
        CellPair pair1(curr_cell,prev_cell);

        // std::cerr << "generated walls between: "
        //           << pair0.first_cell.transpose() << " -> "
        //           << pair0.second_cell.transpose() << std::endl
        //           << pair1.first_cell.transpose() << " -> "
        //           << pair1.second_cell.transpose() << std::endl;
        _map.insert(std::make_pair(pair0, plane_matchable));
        _map.insert(std::make_pair(pair1, plane_matchable));
      }
      ++i;
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

    _trajectory.resize(_num_poses);
    Eigen::Vector3f position(_dimension.x()/2-1+_resolution/2.0,
                             _dimension.y()/2-1+_resolution/2.0,
                             0.0f);

    Eigen::Vector3f traj_point = Eigen::Vector3f::Zero();
    traj_point.z() += 1;

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
      // std::cerr << "finding wall between: " << p.transpose() << " -> " << np.transpose() << std::endl;
      CellPairPlaneMap::iterator it = _map.find(pair);

      if(it != _map.end()){
        // std::cerr << "removed!" << std::endl;
        MatchablePtr m = it->second;
        MatchablePtrSet::iterator jt = _set.find(m);
        if(jt != _set.end())
          _set.erase(jt);
      } else {
        // std::cerr << pair.first_cell.transpose()
                  // << " - " << pair.second_cell.transpose()
                  // << " not found!" << std::endl;
      }

      traj_point.head(2) = position.head(2);
      _trajectory[count] = traj_point;
      position = new_position;
      count++;

      if(count == _num_poses) {
        continue_=false;
      }

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
