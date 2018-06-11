#include "simulator.h"

namespace srrg_g2o_simulator{

  using namespace srrg_matchable;

  void Simulator::populateSet(){
    //ia random number generator
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> uniform_distribution(0, 1);

    //ia insert points
    for (size_t i = 0; i < _num_points; ++i) {
      float x = uniform_distribution(gen)*_dimension.x();
      float y = uniform_distribution(gen)*_dimension.y();
      float z = uniform_distribution(gen);

      MatchablePtr point_matchable(new Matchable(Matchable::Point,
                                                 Eigen::Vector3f(x,y,z)));
      _set.insert(point_matchable);
    }


    //ia new lines generation
    Eigen::Vector2f line_drawing_resolution(_resolution,0.0f);
    for (size_t i = 0; i < _num_lines; ++i) {
      float x = uniform_distribution(gen)*_dimension.x();
      float y = uniform_distribution(gen)*_dimension.y();
      float z = uniform_distribution(gen);

      Eigen::Matrix3f rotation = Eigen::Matrix3f::Identity();
      Eigen::Vector3f axis = Eigen::Vector3f::Zero();


      const double axis_selector = uniform_distribution(gen);
      if (axis_selector < 0.3) {
        axis = Eigen::Vector3f::UnitX();
      } else if (axis_selector < 0.6 && 0.3 < axis_selector) {
        axis = Eigen::Vector3f::UnitY();
      } else {
        axis = Eigen::Vector3f::UnitZ();
      }
      
      MatchablePtr line_matchable(new Matchable(Matchable::Line,
                                                Eigen::Vector3f(x,y,z),
                                                axis,
                                                line_drawing_resolution));
      
      computeRotationMatrix(rotation,axis);
      line_matchable->setRotationMatrix(rotation);
      _set.insert(line_matchable);
    }

    
    // ia new planes generation
    int i = 0;
    uint32_t xcnt = 0;
    uint32_t ycnt = 0;
    uint32_t zcnt = 0;
    while (i < _num_planes) {
      float x = 0;
      float y = 0;
      float z = 0;

      //ia grid coords
      int r = 0, c = 0;
      
      Eigen::Vector3f n=Eigen::Vector3f::Zero();
      Eigen::Vector2i prev_cell = Eigen::Vector2i::Zero();
      Eigen::Vector2i curr_cell = Eigen::Vector2i::Zero();

      const double x_prob = 0.25;
      const double y_prob = 0.25;
      const double z_prob = 1.0 - x_prob - y_prob;

      const double selector = uniform_distribution(gen);
      if (selector < x_prob) {
        //ia x axis
        n = Eigen::Vector3f::UnitX();
        r = round(uniform_distribution(gen)*_dimension.x());
        c = round(uniform_distribution(gen)*_dimension.y());
        
        x = r * _resolution;
        y = c * _resolution + _resolution/2.0;
        
        z = uniform_distribution(gen);
        
        prev_cell.x() = r - 1;
        prev_cell.y() = c;
        curr_cell.x() = r;
        curr_cell.y() = c;
        // std::cerr << "x axis (r,c) = (" << r << ", " << c << ")" << std::endl;
        ++xcnt;
      } else if (selector < x_prob+y_prob && x_prob < selector) {
        //ia y axis
        n = Eigen::Vector3f::UnitY();
        r = round(uniform_distribution(gen)*_dimension.x());
        c = round(uniform_distribution(gen)*_dimension.y());

        x = r * _resolution + _resolution/2.0;
        y = c * _resolution;
          
        z = uniform_distribution(gen);

        prev_cell.x() = r;
        prev_cell.y() = c - 1;
        curr_cell.x() = r;
        curr_cell.y() = c;

        // std::cerr << "y axis (r,c) = (" << r << ", " << c << ")" << std::endl;
        ++ycnt;
      } else {
        // z axis (floor)
        n = Eigen::Vector3f::UnitZ();
        x = round(uniform_distribution(gen)*_dimension.x()) + _resolution/2.0f;
        y = round(uniform_distribution(gen)*_dimension.y()) + _resolution/2.0f;
        ++zcnt;
      }

      // std::cerr << "x: " << xcnt << "\ty: " << ycnt << "\tz: " << zcnt << std::endl;

      if (x - _resolution/2 < 0 || x + _resolution/2 >= _dimension.x() ||
          y - _resolution/2 < 0 || y + _resolution/2 >= _dimension.y() ||
          (z - _resolution/2 < 0 && selector < z_prob)) {
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
      if (selector < z_prob) {
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
      // std::cerr << "***************************************************" << std::endl;
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

    //ia visualization
    _trajectory.resize(_num_poses);
    Eigen::Vector3f traj_point = Eigen::Vector3f::Zero();
    traj_point.z() += 1;

    //ia starting motion
    Eigen::Vector3f pose_2d(_dimension.x()/2-1+_resolution/2.0,
                            _dimension.y()/2-1+_resolution/2.0,
                            0.0f);
    Eigen::Vector3f new_pose_2d = Eigen::Vector3f::Zero();

    while(continue_){

      //sample new position
      float x = 0, y = 0;
      float theta = 0;
      float dir_selector = dis(gen);

      //go forward
      if(dir_selector < 0.5f){
        x = round(cos(pose_2d.z()));
        y = round(sin(pose_2d.z()));
      }
      //go left
      if(dir_selector >= 0.5f && dir_selector < 0.75f){
        x = round(-sin(pose_2d.z()));
        y = round(cos(pose_2d.z()));
        theta = M_PI/2.0f;
      }
      //go right
      if(dir_selector >= 0.75f){
        x = round(sin(pose_2d.z()));
        y = round(-cos(pose_2d.z()));
        theta = -M_PI/2.0f;
      }

      new_pose_2d = pose_2d+Eigen::Vector3f(x,y,theta);
      
      //check if new position is out of grid
      if(new_pose_2d.x() < 0 || new_pose_2d.x() >= _dimension.x()-1 ||
         new_pose_2d.y() < 0 || new_pose_2d.y() >= _dimension.y()-1){
        continue;
      }

      //ia check in the wall map
      const CellPair cell_motion(pose_2d.head(2).cast<int>(), new_pose_2d.head(2).cast<int>());


      //sbraco
      std::cerr << "finding wall between: " << pose_2d.head(2).transpose()
                << " -> " << new_pose_2d.head(2).transpose() << std::endl;
      CellPairPlaneMap::iterator it = _map.find(cell_motion);

      if(it != _map.end()){
        std::cerr << "removed!" << std::endl;
        MatchablePtr m = it->second;
        MatchablePtrSet::iterator jt = _set.find(m);
        if(jt != _set.end())
          _set.erase(jt);
      } else {
        // std::cerr << pair.first_cell.transpose()
                  // << " - " << pair.second_cell.transpose()
                  // << " not found!" << std::endl;
      }

      traj_point.head(2) = pose_2d.head(2);
      _trajectory[count] = traj_point;
      pose_2d = new_pose_2d;
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
