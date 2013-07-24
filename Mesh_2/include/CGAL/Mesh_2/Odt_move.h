
//******************************************************************************
// File Description : Odt move function
//******************************************************************************

#ifndef CGAL_MESH_2_ODT_MOVE_H
#define CGAL_MESH_2_ODT_MOVE_H

#include <string>
#include <set>

namespace CGAL {

namespace Mesh_2 {

template <typename CDT>
class Odt_move
{

  typedef typename CDT::Vertex_handle           Vertex_handle;
  typedef typename CDT::Geom_traits::Vector_2   Vector_2;
  typedef typename CDT::Polygon_2               Polygon_2;
  typedef typename CDT::Point                   Point_2;
  typedef std::set<Point_2>                     Point_set;
  
public:
  /**
   * @brief Return move to apply on \c v according to Odt optimization 
   * function
   */
  Vector_2 operator()(const Vertex_handle& v,
                      CDT& cdt) const
  {
    if (!cdt.cell_is_infinite(v)){
      return CGAL::NULL_VECTOR;
    }
    //to review
    // Calculate the average of \c circumcenters incident to current vertex.
    Polygon_2 poly = cdt.dual(v);
    unsigned int poly_size = poly.size();
    if(poly_size==0)
      return Vector_2(v->point(),v->point());
    else{

      Point_set pset;
      double new_x=0.0;
      double new_y=0.0;

      for(typename Polygon_2::iterator pit=poly.vertices_begin(); pit!=poly.vertices_end(); pit++){
        pset.insert(*pit);
      }
      for(typename Point_set::iterator psit=pset.begin(); psit!=pset.end(); psit++){
        new_x += psit->x();
        new_y += psit->y();
      }

      return Vector_2(v->point(),Point_2(new_x/pset.size(),new_y/pset.size()));
    }
    return CGAL::NULL_VECTOR;
  }
  
#ifdef CGAL_MESH_2_OPTIMIZER_VERBOSE
  static std::string name() { return std::string("Odt"); }
#endif

private:
  Vector_2 lloyd_move_inside_domain(const Vertex_handle& v,
                                    const CDT& cdt) const
  {
    Vector_2 move = CGAL::NULL_VECTOR;
  }
  
};

} // end namespace Mesh_2

} //namespace CGAL

#endif // CGAL_MESH_2_ODT_MOVE_H
