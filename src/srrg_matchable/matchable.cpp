#include "matchable.h"
#include <Eigen/SVD>
#include <Eigen/Eigenvalues>
#include <iostream>

#include "srrg_gl_helpers/opengl_primitives.h"

#include "srrg_image_utils/bresenham.h"


namespace srrg_matchable {

  using namespace srrg_boss;
  using namespace srrg_gl_helpers;
  using namespace srrg_core;
  using namespace std;
  
  Matchable::Matchable(){
  };

  Matchable::Matchable(Matchable::Type type_,
                       const Eigen::Vector3f& point_,
                       const Eigen::Vector3f& direction_vector_,
                       const Eigen::Vector2f &extent_):
    _type(type_),
    _point(point_),
    _direction_vector(direction_vector_),
    _extent(extent_){
    
    _descriptor = 0;
    _direction_matrix.setZero();
    _rotation_matrix.setZero();
    _cum_p.setZero();
    _cum_d.setZero();
    _origin.setZero();
    _cloud.resize(0);
    _contour_cloud.resize(0);
    _association = nullptr;
    _association_id = -1;
    _is_good = true;
    _age = 1;
    _map_id = -1;
    _selected = false;
    updateDirectionMatrix();    
  }

  Matchable* Matchable::clone() const {
    return new Matchable(*this);
  }


  void Matchable::transformInPlace(const Eigen::Isometry3f& T) {
    _point=T*_point;
    if (_type==Point){
      return;
    }
    _direction_vector=T.linear()*_direction_vector;
    _direction_matrix=T.linear()*_direction_matrix*T.linear().transpose();
    //    _rotation_matrix=T.linear()*_rotation_matrix*T.linear().transpose();
    _rotation_matrix=T.linear()*_rotation_matrix;
  }
  
  void Matchable::updateDirectionMatrix(){
    if(_type==Point){
      _direction_matrix.setIdentity();
      return;
    }

    _direction_vector.normalize();

    float d = sqrt(_direction_vector.x()*_direction_vector.x() + _direction_vector.y()*_direction_vector.y());

    const float& dirx = _direction_vector.x();
    const float& diry = _direction_vector.y();
    const float& dirz = _direction_vector.z();

    Eigen::Matrix3f R;
    if(d > std::numeric_limits<float>::min()) {
      R << dirx, diry/d,     dirx*dirz/d,
        diry, -dirx/d,   diry*dirz/d,
        dirz, 0,          -d;
      
      } else {
      R.setIdentity();
    }

    Eigen::DiagonalMatrix<float,3> S;
    S.setZero();
    switch(_type){
    case Line:
      S.diagonal()[1]=1;
      S.diagonal()[2]=1;
      break;
    case Plane:
      S.diagonal()[0]=1;
      break;
    case Surfel:
      S.diagonal()[0]=1;
      S.diagonal()[1]=1e-2f;
      S.diagonal()[2]=1e-2f;
      break;
    default: break;
    }
    _direction_matrix=R*S*R.transpose();
  }

  void Matchable::drawDirectionMatrix() const {
    using namespace srrg_gl_helpers;
    if (_type==Point)
      return;

    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3f> eigensolver;
    eigensolver.compute(_direction_matrix);
    Eigen::Isometry3f transform;
    transform.linear()=eigensolver.eigenvectors();
    transform.translation()=_point;
    float eps=1e-1;
    Eigen::Vector3f inverse_eigenvalues=(eigensolver.eigenvalues()+Eigen::Vector3f(eps, eps, eps)).cwiseInverse();
    inverse_eigenvalues*=0.05;

    glPushMatrix();
    glMultMatrix(transform);
    drawEllipsoid(inverse_eigenvalues(0),
                  inverse_eigenvalues(1),
                  inverse_eigenvalues(2));
    glPopMatrix();
  }

  void Matchable::voxelizeCloud(){
    
    float x_min = std::numeric_limits<float>::max();
    float y_min = std::numeric_limits<float>::max();
    float z_min = std::numeric_limits<float>::max();
    Eigen::Vector3f centroid;
    centroid.setZero();
    for(int i=0; i<_cloud.size(); ++i){
      const Eigen::Vector3f& point = _cloud[i].point();

      if(point.x() < x_min)
        x_min = point.x();
      if(point.y() < y_min)
        y_min = point.y();
      if(point.z() < z_min)
        z_min = point.z();

      centroid += point;
    }

    centroid /= (float)_cloud.size();
    Eigen::Vector3f origin(x_min,y_min,z_min);
    
    float resolution = 0.01f;
    float inverse_resolution = 1.0f/resolution;
    Eigen::Vector3f offset (resolution/2.0f,resolution/2.0f,0.0f);

    std::vector<Eigen::Vector3i, Eigen::aligned_allocator<Eigen::Vector3i> > ipoints(_cloud.size());
    for(int i=0; i<_cloud.size(); ++i){
      const Eigen::Vector3f& point = _cloud[i].point();

      Eigen::Vector3i idx = ((point-(centroid+origin))*inverse_resolution).cast<int>();
      ipoints[i] = idx;
    }

    std::sort(ipoints.begin(),
              ipoints.end(),
              [] (const Eigen::Vector3i& p1, const Eigen::Vector3i& p2) -> bool{
                for (int i=0; i<3; ++i) {
                  if (p1[i]<p2[i])
                    return true;
                  if (p1[i]>p2[i])
                    return false;
                }
                return false;
              });

    Cloud3D* grid_cloud = new Cloud3D;
    grid_cloud->resize(ipoints.size());
    int idx = -1;  
    for (int i=1; i<ipoints.size(); ++i){
      if(idx>=0 && ipoints[i] == ipoints[i-1]){
        continue;
      } else {
        idx++;
        grid_cloud->at(idx) = RichPoint3D((origin+centroid) + ipoints[i].cast<float>()*resolution +offset,Eigen::Vector3f(0,0,0),1.0f);
      }
    }
    grid_cloud->resize(idx);
    const int num_bins=60;
    float angle_res=num_bins/(2*M_PI);
    float angle_ires=1.f/angle_res; 
    std::vector<float> ranges(num_bins);
    std::vector<int> indices(num_bins);
    std::fill(ranges.begin(), ranges.end(), -1.0f);
    std::fill(indices.begin(), indices.end(), -1);

    //simulate a 360 deg scan
    for(int i=0; i < grid_cloud->size(); ++i){
      const Eigen::Vector3f& point = grid_cloud->at(i).point();

      const float range = sqrt(point.x()*point.x() + point.y()*point.y());
      const float angle = atan2(point.x(),point.y())+M_PI;
      int bin=std::round(angle_res*angle);
        
      if (bin<0)
        bin=0;
      if (bin>=num_bins)
        bin=num_bins-1;
      if (ranges[bin]<range) {
        ranges[bin]=range;
        indices[bin]=i;
      }
    }

    delete grid_cloud;

    //fill vertex cloud
    _contour_cloud.clear();
    _contour_cloud.resize(num_bins);
    int k=0;

    for(int i=0; i < num_bins; ++i){

      if(indices[i] == -1) {
        continue;
      }
      float range = ranges[i];
      if(range < 0.05){
        continue;
      }

      float angle = (angle_ires*i)-M_PI;

      Eigen::Vector3f vertex(range*sin(angle), range*cos(angle), 0);
      _contour_cloud[k] = RichPoint3D(vertex,Eigen::Vector3f(0,0,1),1.0f);
      k++;
    }
    _contour_cloud.resize(k);

  }

  void Matchable::update(){

    if(_type == Matchable::Point || _type == Matchable::Surfel){
      _point += _cum_p;
      _point /= (float)_age;
    }

    _direction_vector += _cum_d;
    _direction_vector /= (float)_age;
    _direction_vector.normalize();
    updateDirectionMatrix();
    
  }
  
  void Matchable::merge(const MatchablePtr& m,
                        const Eigen::Isometry3f& T){

    if(_type != m->type())
      throw std::runtime_error("type mismatch!");

    Eigen::Vector3f q;
    float d;
    Eigen::Isometry3f transform,m_transform;

    switch(_type){
    case Matchable::Point:
      _point = (_age*_point + T*m->point())/(float)(_age+1);      
      break;
      
    case Matchable::Line:
      _direction_vector = (_age*_direction_vector + T.linear()*m->directionVector())/(float)(_age+1);      
      _direction_vector.normalize();
      
      q = T*m->point();
      d = (q-_point).dot(_direction_vector);

      m->cloud().transformInPlace(transform.inverse()*m_transform);
      _cloud.add(m->cloud());
      _cloud.voxelize(0.05f);
      
      if(d >= 0){
        if(_extent.x() < d + m->extent().x())
          _extent.x() = d + m->extent().x();
      } else {
        _point = q;
        if(m->extent().x() < _extent.x()-d)
          _extent.x() = _extent.x()-d;
      }
      break;
    case Matchable::Plane:
      _direction_vector = (_age*_direction_vector + T.linear()*m->directionVector())/(float)(_age+1);      
      _direction_vector.normalize();

      transform.linear() = _rotation_matrix;
      transform.translation() = _point;

      m_transform.linear() = T.linear()*m->rotationMatrix();
      m_transform.translation() = T*m->point();
      
      m->cloud().transformInPlace(transform.inverse()*m_transform);

      _cloud.add(m->cloud());
      _cloud.voxelize(0.05f);
      voxelizeCloud();
      break;
    case Matchable::Surfel:
      _point = (_age*_point + T*m->point())/(float)(_age+1);      

      _direction_vector = (_age*_direction_vector + T*m->directionVector())/(float)(_age+1);      
      _direction_vector.normalize();

      break;
    default:
      break;
    }
    _descriptor = m->descriptor();
    _age = std::max(_age,m->age());
    updateDirectionMatrix();    
  }

  void Matchable::getInformation(Eigen::Matrix3f& omega) {
    float mean = 0.f;
    float sigma = 0.f;
    omega.setZero();
    const int cloud_size = _cloud.size();
    if(!cloud_size)
      return;
    
    for(const RichPoint3D& rp : _cloud) {      
      mean += rp.point().norm();
    }
    mean /= (float) cloud_size;
    if(mean > 1e-4f)
      mean = 1.f / mean;
    else
      mean = 1.f;

    const float max = 5.f;
    const float min = 5e-2f;
    const float idiff = 1.f / (mean - min);
    const float pw = (max - min) * idiff;
       
    omega << pw, 0, 0,
      0, pw, 0,
      0, 0, pw;
  }
  
  const void Matchable::weight(Eigen::Matrix3f& w) {
    w.setZero();
    const float max = 5.f;
    const float min = 5e-2f;
    const float d = _point.norm();
    const float idiff = 1.f / (d - min);
    const float pw = (max - min) * idiff;

    switch(_type) {
    case Point:
      w(0,0) = pw;
      w(1,1) = pw;
      w(2,2) = pw;
      return;
      break;
    case Line:
      getInformation(w);
      return;
      break;
    case Plane:
      getInformation(w);
      return;
      break;
    case Surfel:
      w(0,0) = .5f * pw;
      w(1,1) = .5f * pw;
      w(2,2) = .5f * pw;
      return;
      break;
    default:
      getInformation(w);
      return;
      break;
    }
  }

  void Matchable::drawPoint() const {
    if (_type!=Point)
      return;

    Eigen::Isometry3f transform = Eigen::Isometry3f::Identity();
    transform.translation()=_point;

    glPushMatrix();
    glMultMatrix(transform);

    drawSphere(0.05f);
    glPopMatrix();

  }
  
  void Matchable::drawDirectionVector() const {
    if (_type==Point)
      return;
    float line_length= 0.025f;
    glPushAttrib(GL_LINE_WIDTH);
    glLineWidth(2);
    glBegin(GL_LINES);
    glNormal3f(-1.f, -1.f, -1.f);
    glVertex3f(_point.x(), _point.y(), _point.z());
    glVertex3f(_point.x()+_direction_vector.x()*line_length,
               _point.y()+_direction_vector.y()*line_length,
               _point.z()+_direction_vector.z()*line_length);
    glEnd();
    glPopAttrib();
  }

  void Matchable::drawLine() const {
    
    if (_type!=Line)
      return;
    
    float line_length= (_extent.x() > 0.f)? _extent.x(): 10;
    getMatchableColor();

    Eigen::Isometry3f transform;    
    transform.linear()=_rotation_matrix;
    transform.translation()=_point;

    float radius = 0.025;
    float angle;
    float x,y;
    float angle_stepsize = 0.017;
    
    glPushMatrix();
    glMultMatrix(transform);

    glBegin(GL_QUAD_STRIP);
    angle = 0.0;
    while( angle < 2*M_PI ) {
      x = radius * cos(angle);
      y = radius * sin(angle);
      glVertex3f(x,y,line_length);
      glVertex3f(x,y,0.0);
      angle = angle + angle_stepsize;
    }
    glVertex3f(radius, 0.0, line_length);
    glVertex3f(radius, 0.0, 0.0);
    glEnd();
    
    glPopMatrix();

  }

  void Matchable::drawPlane() const {

    if (_type!=Plane)
      return;

    Eigen::Isometry3f transform;    
    transform.linear()=_rotation_matrix;
    transform.translation()=_point;
    
    glPushMatrix();
    glMultMatrix(transform);
    getMatchableColor();
//    glBegin(GL_LINE_LOOP);
//    glBegin(GL_TRIANGLE_FAN);

//    glVertex3f(0,0,0);
//    glNormal3f(0,0,-1);

//    for(int i=_contour_cloud.size()-1; i > 0 ; --i){
//    glVertex3f(_contour_cloud[i].point().x(),
//               _contour_cloud[i].point().y(),
//               _contour_cloud[i].point().z());
//    glNormal3f(0,0,-1);
//    }

//    glVertex3f(_contour_cloud[_contour_cloud.size()-1].point().x(),
//               _contour_cloud[_contour_cloud.size()-1].point().y(),
//               _contour_cloud[_contour_cloud.size()-1].point().z());
//    glNormal3f(0,0,-1);

//    glEnd();

    glBegin(GL_TRIANGLE_FAN);

    glVertex3f(0,0,0);
    glNormal3f(0,0,1);

    glVertex3f(_extent.x(),_extent.y(),0);
    glNormal3f(0,0,1);

    glVertex3f(-_extent.x(),_extent.y(),0);
    glNormal3f(0,0,1);

    glVertex3f(-_extent.x(),-_extent.y(),0);
    glNormal3f(0,0,1);

    glVertex3f(_extent.x(),-_extent.y(),0);
    glNormal3f(0,0,1);

    glVertex3f(_extent.x(),_extent.y(),0);
    glNormal3f(0,0,1);

    glEnd();

    glPopMatrix();

  }


 
  void Matchable::draw(bool show_direction_matrix,
                       const Eigen::Vector3f& base_color) const{

    glPushAttrib(GL_COLOR|GL_POINT_SIZE|GL_LINE_WIDTH);

    getMatchableColor();

    if(_selected)
      glPointSize(10);

    drawPoint();

    drawDirectionVector();

    drawLine();
    drawPlane();
    if (show_direction_matrix)
      drawDirectionMatrix();
    glPopAttrib();
  }

  void Matchable::getMatchableColor() const {
    if(_map_id == -1)
      switch(_type){
        case Matchable::Point: glColor4f(1.0f,0.0f,1.0f,1.f); break;
        case Matchable::Line: glColor4f(1.0f,1.0f,0.0f,1.f); break;
        case Matchable::Plane: glColor4f(0.0f,1.0f,1.0f,1.f); break;
      case Matchable::Surfel: glColor4f(1.0f,1.0f,0.0f,1.f); break;
      }
    else
      switch(_type){
      case Matchable::Point: glColor4f(1.0f,0.0f,0.0f,1.0f); break;
      case Matchable::Line: glColor4f(0.0f,1.0f,0.0f,0.75f); break;
      case Matchable::Plane: glColor4f(0.0f,0.0f,1.0f,.95f); break;
      case Matchable::Surfel: glColor4f(1.0f,0.0f,1.0f,0.5f); break;

      }      
  }

  void Matchable::serialize(ObjectData& data_, IdContext& context) {
    data_.setInt("type", _type);
    SRRG_TO_BOSS_MATRIX(data_, point);
    SRRG_TO_BOSS_MATRIX(data_, direction_vector);
    SRRG_TO_BOSS_MATRIX(data_, direction_matrix);
    SRRG_TO_BOSS_MATRIX(data_, extent);
    SRRG_TO_BOSS_MATRIX(data_, rotation_matrix);
    data_.setInt("age", _age);
    Eigen::MatrixXi _descriptor_matrix;
    cv2eigen(_descriptor, _descriptor_matrix);
    SRRG_TO_BOSS_MATRIX(data_, descriptor_matrix);

  }

  void Matchable::deserialize(ObjectData& data_, IdContext& context) {
    int type = data_.getInt("type");
    switch(type) {
    case 0:
      _type = Point;
      break;
    case 1:
      _type = Line;
      break;
    case 2:
      _type = Plane;
      break;
    case 3:
      _type = Surfel;
      break;
    default:
      break;
    }
    SRRG_FROM_BOSS_MATRIX(data_, point);
    SRRG_FROM_BOSS_MATRIX(data_, direction_vector);
    SRRG_FROM_BOSS_MATRIX(data_, direction_matrix);
    SRRG_FROM_BOSS_MATRIX(data_, extent);
    SRRG_FROM_BOSS_MATRIX(data_, rotation_matrix);
    _age = data_.getInt("age");

    Eigen::MatrixXi _descriptor_matrix;
    SRRG_TO_BOSS_MATRIX(data_, descriptor_matrix);
    eigen2cv(_descriptor_matrix, _descriptor);

  }
}// end namespace


