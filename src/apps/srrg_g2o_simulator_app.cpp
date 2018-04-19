#include <iostream>
#include <srrg_g2o_simulator_viewer/simulator_viewer.h>

#include <qapplication.h>
#include <qglviewer.h>

using namespace srrg_core;
using namespace srrg_matchable;
using namespace srrg_g2o_simulator;


void computeRotationMatrix(Eigen::Matrix3f& rotation_matrix,
                           const Eigen::Vector3f& direction){

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

bool neighborVisited(const std::vector<std::vector<int> > &visited,const Eigen::Vector2i &p, int x_dim, int y_dim){
  const Eigen::Vector2i p_right (p.x()+1,p.y());
  if(!(p_right.x() < 0 || p_right.x() >= x_dim-1 ||
     p_right.y() < 0 || p_right.y() >= y_dim-1)){
    if(!visited[p_right.x()][p_right.y()])
      return false;
  }
  const Eigen::Vector2i p_up (p.x(),p.y()+1);
  if(!(p_up.x() < 0 || p_up.x() >= x_dim-1 ||
     p_up.y() < 0 || p_up.y() >= y_dim-1)){
    if(!visited[p_up.x()][p_up.y()])
      return false;
  }
  const Eigen::Vector2i p_left (p.x()-1,p.y());
  if(!(p_left.x() < 0 || p_left.x() >= x_dim-1 ||
     p_left.y() < 0 || p_left.y() >= y_dim-1)){
    if(!visited[p_left.x()][p_left.y()])
      return false;
  }
  const Eigen::Vector2i p_down (p.x(),p.y()-1);
  if(!(p_down.x() < 0 || p_down.x() >= x_dim-1 ||
     p_down.y() < 0 || p_down.y() >= y_dim-1)){
    if(!visited[p_down.x()][p_down.y()])
      return false;
  }
  return true;
}

struct CellPair{
  CellPair(const Eigen::Vector2i &first_cell_,
           const Eigen::Vector2i &second_cell_):
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
typedef std::map<CellPair,MatchablePtr> CellPairPlaneMap;
typedef std::set<MatchablePtr> MatchablePtrSet;

int main(int argc, char **argv){

  float resolution = 1.0f;
  Eigen::Vector2i dimension (10,10);

  MatchablePtrSet set;
  CellPairPlaneMap map;

  //insert points
  for(size_t j=0; j<dimension.y(); ++j)
    for(size_t i=0; i<dimension.x(); ++i){
      const Eigen::Vector3f p(i*resolution, j*resolution, 0.0f);
      MatchablePtr point_matchable(new Matchable(Matchable::Point,p));
      set.insert(point_matchable);
    }

  //insert lines along x axes
  for(size_t j=0; j<dimension.y(); ++j) //
    for(size_t i=0; i<dimension.x()-1; ++i){
      const Eigen::Vector3f p(i*resolution, j*resolution, 0.0f);
      const Eigen::Vector3f n(1.0f,0.0f,0.0f);
      MatchablePtr line_matchable(new Matchable(Matchable::Line,p,n,Eigen::Vector2f(resolution,0.0f)));
      Eigen::Matrix3f R;
      computeRotationMatrix(R,n);
      line_matchable->setRotationMatrix(R);
      set.insert(line_matchable);
    }

  //insert lines along y axes
  for(size_t i=0; i<dimension.x(); ++i)
    for(size_t j=0; j<dimension.y()-1; ++j){
      const Eigen::Vector3f p(i*resolution, j*resolution, 0.0f);
      const Eigen::Vector3f n(0.0f,1.0f,0.0f);
      MatchablePtr line_matchable(new Matchable(Matchable::Line,p,n,Eigen::Vector2f(resolution,0.0f)));
      Eigen::Matrix3f R;
      computeRotationMatrix(R,n);
      line_matchable->setRotationMatrix(R);
      set.insert(line_matchable);
    }

  //insert lines along z axes
  for(size_t j=0; j<dimension.y(); ++j)
    for(size_t i=0; i<dimension.x(); ++i){
      const Eigen::Vector3f p(i*resolution, j*resolution, 0.0f);
      const Eigen::Vector3f n(0.0f,0.0f,1.0f);
      MatchablePtr line_matchable(new Matchable(Matchable::Line,p,n,Eigen::Vector2f(resolution,0.0f)));
      Eigen::Matrix3f R;
      computeRotationMatrix(R,n);
      line_matchable->setRotationMatrix(R);
      set.insert(line_matchable);
    }

  //insert planes along x axes
  for(size_t j=0; j<dimension.y()-1; ++j)
    for(size_t i=0; i<dimension.x(); ++i){
      const Eigen::Vector2i cell(i,j);
      const Eigen::Vector3f p(cell.x()*resolution, cell.y()*resolution+resolution/2.0f, resolution/2.0f);
      const Eigen::Vector3f n(1.0f,0.0f,0.0f);
      MatchablePtr plane_matchable(new Matchable(Matchable::Plane,p,n,Eigen::Vector2f(resolution/2.0f,resolution/2.0f)));
      Eigen::Matrix3f R;
      computeRotationMatrix(R,n);
      plane_matchable->setRotationMatrix(R);
      set.insert(plane_matchable);

      if(i >= 0 && i < dimension.x()-1){
        const Eigen::Vector2i left_cell(i-1,j);
        CellPair left_pair(left_cell,cell);
        map[left_pair] = plane_matchable;
        CellPair right_pair(cell,left_cell);
        map[right_pair] = plane_matchable;
      }
    }

  //insert planes along y axes
  for(size_t i=0; i<dimension.x()-1; ++i)
    for(size_t j=0; j<dimension.y(); ++j){
      const Eigen::Vector2i cell(i,j);
      const Eigen::Vector3f p(cell.x()*resolution+resolution/2.0f, cell.y()*resolution, resolution/2.0f);
      const Eigen::Vector3f n(0.0f,1.0f,0.0f);
      MatchablePtr plane_matchable(new Matchable(Matchable::Plane,p,n,Eigen::Vector2f(resolution/2.0f,resolution/2.0f)));
      Eigen::Matrix3f R;
      computeRotationMatrix(R,n);
      plane_matchable->setRotationMatrix(R);
      set.insert(plane_matchable);

      if(j >= 0 && j < dimension.y()-1){
        const Eigen::Vector2i down_cell(i,j-1);
        CellPair down_pair(down_cell,cell);
        map[down_pair] = plane_matchable;
        CellPair up_pair(cell,down_cell);
        map[up_pair] = plane_matchable;
      }
    }

  //insert planes along z axes
  for(size_t j=0; j<dimension.y()-1; ++j)
    for(size_t i=0; i<dimension.x()-1; ++i){
      const Eigen::Vector3f p(i*resolution+resolution/2.0f, j*resolution+resolution/2.0f, 0.0f);
      const Eigen::Vector3f n(0.0f,0.0f,1.0f);
      MatchablePtr plane_matchable(new Matchable(Matchable::Plane,p,n,Eigen::Vector2f(resolution/2.0f,resolution/2.0f)));
      Eigen::Matrix3f R;
      computeRotationMatrix(R,n);
      plane_matchable->setRotationMatrix(R);
      set.insert(plane_matchable);

    }

  std::cerr << "Scene has " << set.size() << " matchables before sbraco" << std::endl;

  //sbraco
  std::vector<std::vector<int> > visited;
  visited.resize(dimension.x());
  for(int i=0; i<dimension.x(); ++i){
    visited[i].resize(dimension.y());
    for(int j=0; j<dimension.y(); ++j)
      visited[i][j] = 0;
  }

  std::cerr << "sbraco..." << std::endl;

  Eigen::Vector3f position(4.0,4.0,M_PI/2.0f);
  visited[4][4]=1;
  const int num_poses=50;
  int count=0;

  bool continue_=true;
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<> dis(0, 1);

  while(continue_){

    std::cerr << std::endl << "current position: " << position.transpose() << std::endl;

    //sample new position
    Eigen::Vector3f increment = Eigen::Vector3f::Zero();
    float n = dis(gen);

    //go forward
    if(n < 0.5f){
      increment.x() = cos(position.z());
      increment.y() = sin(position.z());
      std::cerr << "go forward" << std::endl;
    }
    //go left
    if(n >= 0.5f && n < 0.75f){
      increment.x() = -sin(position.z());
      increment.y() = cos(position.z());
      increment.z() = M_PI/2.0f;
      std::cerr << "go left" << std::endl;
    }
    //go right
    if(n >= 0.75f){
      increment.x() = sin(position.z());
      increment.y() = -cos(position.z());
      increment.z() = -M_PI/2.0f;
      std::cerr << "go right" << std::endl;
    }

    Eigen::Vector3f new_position = position+increment;
    std::cerr << "new position: " << new_position.transpose() << std::endl;

    //check if new position is out of grid
    if(new_position.x() < 0.0f || new_position.x() >= (float)(dimension.x()-1) ||
       new_position.y() < 0.0f || new_position.y() >= (float)(dimension.y()-1)){
      std::cerr << "out of grid!" << std::endl;
      continue;
    }

    //check if new position is already visited
    if(!visited[new_position.x()][new_position.y()]){

      visited[new_position.x()][new_position.y()] = 1;

      std::cerr << "eddaje!" << std::endl;
      const Eigen::Vector2i p = position.head(2).cast<int>();
      const Eigen::Vector2i np = new_position.head(2).cast<int>();
      CellPair pair(p,np);
      CellPairPlaneMap::iterator it = map.find(pair);

      if(it != map.end()){
        MatchablePtr m = it->second;
        MatchablePtrSet::iterator jt = set.find(m);
        if(jt != set.end())
          set.erase(jt);
      } else {
        std::cerr << pair.first_cell.transpose() << " - " << pair.second_cell.transpose() << " not found!" << std::endl;
//        continue;
      }

    }

//    if(neighborVisited(visited,new_position.head(2).cast<int>(),dimension.x(),dimension.y())){
//      std::cerr << "deadlock" << std::endl;
//      continue_=false;
//    }

    position = new_position;
    count++;

    if(count == num_poses)
      continue_=false;

  }

  Scene* scene = new Scene();
  for(MatchablePtrSet::iterator it = set.begin(); it != set.end(); ++it){
    MatchablePtr m = *it;
    scene->push_back(m);
  }

  std::cerr << "Scene has " << scene->size() << " matchables after sbraco" << std::endl;

  QApplication app(argc, argv);
  SimulatorViewer viewer;
  viewer.setScene(scene);
  viewer.show();
  app.exec();

  return 0;
}
