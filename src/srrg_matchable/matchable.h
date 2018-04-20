#pragma once 
#include <srrg_boss/serializable.h>
#include <srrg_boss/eigen_boss_plugin.h>

#include <srrg_types/defs.h>
#include <srrg_types/cloud_3d.h>

// cv2Eigen
#include <opencv2/core/eigen.hpp> 

namespace srrg_matchable {

  class Matchable;
  
  //bdc,  Welcome C++11, damn!
  typedef std::shared_ptr<Matchable> MatchablePtr;
  typedef std::set<MatchablePtr> MatchablePtrSet;

  class Matchable : public srrg_boss::Serializable {
  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
    
    enum Type {Point=0, Line=1, Plane=2, Surfel=3};
    static const int Dim = 4; //< update Dim consistently with Type number
    
    Matchable(Type type,
              const Eigen::Vector3f& point_,
              const Eigen::Vector3f& direction_vector_=Eigen::Vector3f(0,0,0),
              const Eigen::Vector2f& extent_ = Eigen::Vector2f::Zero());

    // Matchable(const Matchable& m);// copy constructor
    
    inline Type type() const {return _type;}
    inline const Eigen::Vector3f& directionVector() const {return _direction_vector;}
    inline const Eigen::Matrix3f& directionMatrix() const {return _direction_matrix;}
    inline const Eigen::Vector3f& point() const {return _point;}
    inline const Eigen::Vector2f& extent() const {return _extent;}
    inline const Eigen::Vector3f& origin() const {return _origin;}
    inline const srrg_core::Cloud3D& cloud() const {return _cloud;}
    inline srrg_core::Cloud3D& cloud() {return _cloud;}
    inline const cv::Mat& descriptor() const {return _descriptor;}
    const inline Eigen::Matrix3f& rotationMatrix() const {return _rotation_matrix;}
    inline Eigen::Vector3f& cumP() {return _cum_p;}
    inline Eigen::Vector3f& cumD() {return _cum_d;}

    const void weight(Eigen::Matrix3f&);
    void getInformation(Eigen::Matrix3f&);
    
    void update();
    
    inline void setPoint(const Eigen::Vector3f& point_){_point = point_;}
    inline void setDirectionVector(const Eigen::Vector3f& direction_vector_){_direction_vector = direction_vector_;}
    inline void setCloud(const srrg_core::Cloud3D& cloud_){_cloud = cloud_;}
    inline void setExtent(const Eigen::Vector2f& extent_){_extent = extent_;}
    inline void setOrigin(const Eigen::Vector3f& origin_){_origin = origin_;}
    inline void setDescriptor(const cv::Mat& descriptor_){_descriptor = descriptor_;}
    inline void setRotationMatrix(const Eigen::Matrix3f& rotation_matrix_){_rotation_matrix = rotation_matrix_;}
    
    void transformInPlace(const Eigen::Isometry3f& T);

    void merge(const MatchablePtr& m,
               const Eigen::Isometry3f& T = Eigen::Isometry3f::Identity());

    //! the one below requires gl. It is exceptionally inefficient :)
    void draw(bool show_direction_matrix=false,
              const Eigen::Vector3f& base_color = Eigen::Vector3f::Zero()) const;

    //! returns a (proper) copy of itself
    Matchable* clone() const;

    //! returns a pointer to the associated (if any) matchable
    inline void setAssociation(MatchablePtr association_,
                               const int association_id_ = -1) {
      _association = association_;
      _association_id = association_id_;
    }
    MatchablePtr association() {return _association;}
    const MatchablePtr association() const {return _association;}
    
    int associationID() {return _association_id;}
    const int associationID() const {return _association_id;}
    
    const bool isGood() const {return _is_good;}
    void setIsGood(bool is_good_){_is_good = is_good_;}
    
    int dof() const {return _type;}
    
    void updateDirectionMatrix();
    const int age() const {return _age;}
    int& age(){return _age;}

    const int mapId() const {return _map_id;}
    int& mapId(){return _map_id;}

    void getMatchableColor() const;

    virtual void serialize(srrg_boss::ObjectData& data, srrg_boss::IdContext& context);
    virtual void deserialize(srrg_boss::ObjectData& data, srrg_boss::IdContext& context);

    void voxelizeCloud();

    inline void select() {_selected = true;}
    inline void unSelect() {_selected = false;}
    
  protected:
    void drawDirectionVector() const;
    void drawDirectionMatrix() const;
    void drawPoint() const;
    void drawLine() const;
    void drawPlane() const;
    Matchable();
  private:
    Type _type;
    cv::Mat _descriptor;
    Eigen::Vector3f _point;
    Eigen::Vector3f _direction_vector;
    Eigen::Matrix3f _direction_matrix;
    Eigen::Vector2f _extent;
    Eigen::Vector3f _origin;
    srrg_core::Cloud3D _cloud;
    srrg_core::Cloud3D _contour_cloud;
    Eigen::Matrix3f _rotation_matrix;
    Eigen::Vector3f _cum_p;
    Eigen::Vector3f _cum_d;
    MatchablePtr _association;
    int _association_id;
    bool _is_good;
    int _age;
    int _map_id;
    bool _selected;
    // if you add something here
    // UPLOAD THE COPY CONSTRUCTOR
  };

} // end namespace srrg_matchable
