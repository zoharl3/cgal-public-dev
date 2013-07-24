
//******************************************************************************
// File Description : Lloyd move function
//******************************************************************************

#ifndef CGAL_MESH_2_LLOYD_MOVE_H
#define CGAL_MESH_2_LLOYD_MOVE_H

#include <string>
#include <CGAL/centroid.h>

namespace CGAL {

namespace Mesh_2 {

template <typename CDT>
class Lloyd_move
{
  typedef typename CDT::Vertex_handle         Vertex_handle;
  typedef typename CDT::Geom_traits::Vector_2 Vector_2;
  typedef typename CDT::Polygon_2             Polygon_2;

public:
  /**
   * @brief Return move to apply on \c v according to Lloyd optimization 
   * function
   */
  Vector_2 operator()(const Vertex_handle& v,
                       CDT& cdt) const
  {
    if (!cdt.cell_is_infinite(v)){
      return CGAL::NULL_VECTOR;
    }
    //to review
    //Vector_2 move = CGAL::NULL_VECTOR;
    Polygon_2 poly = cdt.dual(v);

    // If the dual returns no points,\c that vertex doesn't have dual.
    if(poly.size()==0)
      return Vector_2(v->point(),v->point());

    // Return the vector_2 of the origin\c and the centroid of its voronoi cell.
    return Vector_2(v->point(),CGAL::centroid(poly.vertices_begin(),poly.vertices_end(),CGAL::Dimension_tag<0>()));
  }
  
#ifdef CGAL_MESH_2_OPTIMIZER_VERBOSE
  static std::string name() { return std::string("Lloyd"); }
#endif
  
};

} // end namespace Mesh_2

} //namespace CGAL

#endif // CGAL_MESH_2_LLOYD_MOVE_H
