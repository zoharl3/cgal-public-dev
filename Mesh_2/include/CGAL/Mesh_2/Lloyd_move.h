
//******************************************************************************
// File Description : Lloyd move function
//******************************************************************************

#ifndef CGAL_MESH_2_LLOYD_MOVE_H
#define CGAL_MESH_2_LLOYD_MOVE_H

#include <string>
#include <CGAL/centroid.h>
#include <CGAL/Polygon_2_algorithms.h>
#include <CGAL/intersections.h>

namespace CGAL {

namespace Mesh_2 {

template <typename CDT>
class Lloyd_move
{
  typedef typename CDT::Geom_traits           Gt;
  typedef typename Gt::Segment_2              Segment_2;
  typedef typename CDT::Vertex_handle         Vertex_handle;
  typedef typename CDT::Geom_traits::Vector_2 Vector_2;
  typedef typename CDT::Polygon_2             Polygon_2;
  typedef typename CDT::Point                 Point_2;
  typedef typename CDT::Edge_circulator       Edge_circulator;
  typedef typename CDT::Edge                  Edge;
  typedef typename CDT::Face_handle           Face_handle;
  typedef typename CDT::Face_circulator       Face_circulator;

public:
  /**
   * @brief Return move to apply on \c v according to Lloyd optimization 
   * function
   */
  Vector_2 operator()(const Vertex_handle& v,
                       CDT& cdt) const
  {
    if (cdt.cell_is_infinite(v))
      return CGAL::NULL_VECTOR;
    
    if(v->input_constraint())
      return CGAL::NULL_VECTOR;
    // CHECK IF THERE ARE ONLY TWO CONSTRAINTS WITH THAT VERTEX AND
    // IF YES, THEN CHECK IF THOSE CONSTRAINTS ARE PARALLELS.

    //to review
    //Vector_2 move = CGAL::NULL_VECTOR;
    Polygon_2 poly = cdt.dual(v);
    unsigned int poly_size = poly.size();

    // If the dual returns no points,\c that vertex doesn't have dual.
    if(poly_size==0)
      return CGAL::NULL_VECTOR;

    // Return the vector_2 of the origin\c and the centroid of its voronoi cell.
    //std::cout<<"NEW POSITION"<<std::endl;
    Point_2 new_point = CGAL::centroid(poly.vertices_begin(),
      poly.vertices_end(),
      CGAL::Dimension_tag<0>());
    //std::cout<<poly<<std::endl;
    //std::cout<<"new_point: "<<new_point<<" -> "<<cdt.istrue();

    Edge parallel_constraint;

    Edge_circulator ec=cdt.incident_edges(v), done(ec);
    unsigned int num_constraints = 0;
    do{
      if(cdt.is_constrained(*ec)){
        parallel_constraint = *ec;
        num_constraints++;
      }
      if(num_constraints>2){
        return CGAL::NULL_VECTOR;
      }
      ec++;
    }while(ec != done);

    Segment_2 new_segment = Segment_2(new_point,v->point());

    unsigned int moving = 0;
    Face_handle fh = cdt.locate(new_point);

    switch(num_constraints){
      case 0: moving = 2; break;
      case 1: moving = 0; break;
      case 2: CGAL::parallel(cdt.segment(parallel_constraint),new_segment)?moving=1:moving=0; break;
      default: moving = 0; break;
    }
    if(moving == 2){
      bool moving2 = true;

      Face_circulator face = cdt.incident_faces(v);
      Face_circulator begin = face;
      CGAL_For_all(face, begin){
        if(cdt.is_infinite(face))
          return CGAL::NULL_VECTOR;
        else{
          int vertex_index;
          for(int i=0;i<3;i++)
            if(face->vertex(i) == v){
              vertex_index = i;
              break;
            }

          Edge edge = Edge(face,vertex_index);
          if (cdt.is_constrained(edge)){
            if(CGAL::do_intersect(new_segment,cdt.segment(edge)))
              return CGAL::NULL_VECTOR;
          }
        }
      }
      return Vector_2(v->point(),new_point);  
    }
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
