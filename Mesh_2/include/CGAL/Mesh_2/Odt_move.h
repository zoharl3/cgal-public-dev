
//******************************************************************************
// File Description : Odt move function
//******************************************************************************

#ifndef CGAL_MESH_2_ODT_MOVE_H
#define CGAL_MESH_2_ODT_MOVE_H

#include <string>
#include <set>
#include <CGAL/Polygon_2_algorithms.h>

namespace CGAL {

namespace Mesh_2 {

template <typename CDT>
class Odt_move
{
  typedef typename CDT::Geom_traits             Gt;
  typedef typename Gt::Segment_2                Segment_2;
  typedef typename CDT::Vertex_handle           Vertex_handle;
  typedef typename CDT::Geom_traits::Vector_2   Vector_2;
  typedef typename CDT::Polygon_2               Polygon_2;
  typedef typename CDT::Point                   Point_2;
  typedef std::set<Point_2>                     Point_set;
  typedef typename CDT::Edge_circulator         Edge_circulator;
  typedef typename CDT::Edge                    Edge;
  typedef typename CDT::Face_handle             Face_handle;
  
public:
  /**
   * @brief Return move to apply on \c v according to Odt optimization 
   * function
   */
  Vector_2 operator()(const Vertex_handle& v,
                      CDT& cdt) const
  {
    if (cdt.cell_is_infinite(v)){
      return CGAL::NULL_VECTOR;
    }
    if(v->input_constraint())
      return CGAL::NULL_VECTOR;
    //to review
    // Calculate the average of \c circumcenters incident to current vertex.
    Polygon_2 poly = cdt.dual(v);
    unsigned int poly_size = poly.size();
    if(poly_size==0)
      return CGAL::NULL_VECTOR;
  
    //if(CGAL::bounded_side_2(poly.vertices_begin(),poly.vertices_end(),))

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

    Point_2 new_point = Point_2(new_x/pset.size(),new_y/pset.size());

    Edge parallel_constraint;

    Edge_circulator ec=cdt.incident_edges(v), done(ec);
    unsigned int num_constraints = 0;
    do{
      if(cdt.is_constrained(*ec)){
        parallel_constraint = *ec;
        num_constraints++;
      }
      if(num_constraints>2){// if more than two constraints, it SHOULDN'T MOVE
        return CGAL::NULL_VECTOR;
      }
      ec++;
    }while(ec != done);

    Segment_2 new_segment = Segment_2(new_point,v->point());
    switch(num_constraints){
      case 0: 
        if(cdt.is_inside_triangulation_cell(v,new_point)){
          typename CDT::Locate_type loc;
          int li;
          Face_handle fh = cdt.locate(new_point, loc, li);
          if( loc == CDT::VERTEX ){
            // Delete vertex from here, or just return NULL?
            return Vector_2(v->point(),new_point);  
          }else{
            if(!fh->blind())
              return Vector_2(v->point(),new_point);
            else
              return CGAL::NULL_VECTOR;
          }
        }
        else
          return CGAL::NULL_VECTOR;

      case 1: return CGAL::NULL_VECTOR;
      case 2: 
        if( CGAL::parallel(cdt.segment(parallel_constraint),new_segment)){
          Face_handle fh = cdt.locate(new_point);
          return Vector_2(v->point(),new_point);
        }
        else
          return CGAL::NULL_VECTOR;
      default: return CGAL::NULL_VECTOR;
    }
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
