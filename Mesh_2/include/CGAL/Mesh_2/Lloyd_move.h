
//******************************************************************************
// File Description : Lloyd move function
//******************************************************************************

#ifndef CGAL_MESH_2_LLOYD_MOVE_H
#define CGAL_MESH_2_LLOYD_MOVE_H

#include <string>

namespace CGAL {

namespace Mesh_2 {

template <typename CDT>
class Lloyd_move
{
  typedef typename CDT::Vertex_handle Vertex_handle;
  typedef typename CDT::Geom_traits::Vector_3 Vector_3;

  /**
   * @brief Return move to apply on \c v according to Lloyd optimization 
   * function
   */
  Vector_3 operator()(const Vertex_handle& v,
                      const CDT& cdt) const
  {
    //todo
    return CGAL::NULL_VECTOR;
  }
  
#ifdef CGAL_MESH_2_OPTIMIZER_VERBOSE
  static std::string name() { return std::string("Lloyd"); }
#endif
  
};

} // end namespace Mesh_2

} //namespace CGAL

#endif // CGAL_MESH_2_LLOYD_MOVE_H
