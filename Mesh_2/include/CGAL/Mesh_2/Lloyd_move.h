
//******************************************************************************
// File Description : Lloyd move function
//******************************************************************************

#ifndef CGAL_MESH_2_LLOYD_MOVE_H
#define CGAL_MESH_2_LLOYD_MOVE_H

#include <string>
#include <CGAL/centroid.h>
#include <CGAL/Polygon_2_algorithms.h>

namespace CGAL {

namespace Mesh_2 {

template <typename CDT>
class Lloyd_move
{
  typedef typename CDT::Geom_traits           Gt;
  typedef typename CDT::Vertex_handle         Vertex_handle;
  typedef typename CDT::Geom_traits::Vector_2 Vector_2;
  typedef typename CDT::Polygon_2             Polygon_2;
  typedef typename CDT::Point                 Point_2;

public:
  /**
   * @brief Return move to apply on \c v according to Lloyd optimization 
   * function
   */
  Vector_2 operator()(const Vertex_handle& v,
                       CDT& cdt) const
  {
    if (cdt.cell_is_infinite(v)){
      //std::cout<<"INFINITE"<<std::endl;
      return CGAL::NULL_VECTOR;
    }
    //to review
    //Vector_2 move = CGAL::NULL_VECTOR;
    Polygon_2 poly = cdt.dual(v);

    // If the dual returns no points,\c that vertex doesn't have dual.
    if(poly.size()==0){
      //std::cout<<"SIZE 0"<<std::endl;
      return CGAL::NULL_VECTOR;
    }
    // Return the vector_2 of the origin\c and the centroid of its voronoi cell.
    //std::cout<<"NEW POSITION"<<std::endl;
    Point_2 new_point = CGAL::centroid(poly.vertices_begin(),poly.vertices_end(),CGAL::Dimension_tag<0>());
    //std::cout<<poly<<std::endl;
    //std::cout<<"new_point: "<<new_point<<" -> "<<cdt.istrue();

    if(cdt.is_inside_triangulation_cell(v,new_point))
      return Vector_2(v->point(),new_point);
    else
      return CGAL::NULL_VECTOR;

    //Point_2 points[] = { Point_2(-4,0), Point_2(0,-1), Point_2(4,0), Point_2(2,0.6), Point_2(0,1)};
    /*if(CGAL::bounded_side_2(points,points+5,new_point, Gt())== CGAL::ON_BOUNDED_SIDE) {
      return Vector_2(v->point(),new_point);
    }else{
      return CGAL::NULL_VECTOR;
    }*/

    /*if(CGAL::bounded_side_2(points,points+5,new_point, Gt())) {
      case CGAL::ON_BOUNDED_SIDE :
        std::cout << " is inside the polygon.\n";
        return Vector_2(v->point(),new_point);
        break;
      case CGAL::ON_BOUNDARY:
        std::cout << " is on the polygon boundary.\n";
        break;
      case CGAL::ON_UNBOUNDED_SIDE:
        std::cout << " is outside the polygon.\n";
        break;
    }*/
    /*if(CGAL::bounded_side_2(poly.vertices_begin(),poly.vertices_end(),new_point,Gt())
      == CGAL::ON_BOUNDED_SIDE){*/
      //return Vector_2(v->point(),new_point);
    /*}else{
      return CGAL::NULL_VECTOR;
    }*/
  }
  
#ifdef CGAL_MESH_2_OPTIMIZER_VERBOSE
  static std::string name() { return std::string("Lloyd"); }
#endif
  
};

} // end namespace Mesh_2

} //namespace CGAL

#endif // CGAL_MESH_2_LLOYD_MOVE_H
