// Copyright (c) 2009 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
//
//
// Author(s)     : Raul Gallegos, Jane Tournois
//
//******************************************************************************
// File Description : Odt move function
//******************************************************************************

#ifndef CGAL_MESH_2_ODT_MOVE_H
#define CGAL_MESH_2_ODT_MOVE_H

#include <string>
#include <set>
#include <CGAL/Polygon_2_algorithms.h>
#include <CGAL/Mesh_3/Uniform_sizing_field.h>

namespace CGAL {

namespace Mesh_2 {

template <typename CDT,
  typename SizingField = Uniform_sizing_field<typename CDT::Triangulation> >
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
  typedef SizingField Sizing_field;
  /**
   * @brief Return move to apply on \c v according to Odt optimization 
   * function
   */
  Vector_2 operator()(const Vertex_handle& v,
                      CDT& cdt,
                       const Sizing_field& sizing_field = Sizing_field() ) const
  {
    //TODO: check if infinite vertices should be treated like this.
    if (cdt.cell_is_infinite(v)){
      return CGAL::NULL_VECTOR;
    }

    // If a vertex is part of and input constraint, then it shouldn't move.
    if(v->input_constraint())
      return CGAL::NULL_VECTOR;

    // Calculate the average of circumcenters incident to current vertex.
    Polygon_2 poly = cdt.dual(v);
    std::size_t poly_size = poly.size();

    CGAL_assertion(poly_size!=0);

    double new_x=0.0;
    double new_y=0.0;

    // Getting average point from the dual
    for(typename Polygon_2::Vertex_iterator psit=poly.vertices_begin(); psit!=poly.vertices_end(); psit++){
      new_x += psit->x();
      new_y += psit->y();
    }

    Point_2 new_point = Point_2(new_x/poly_size,new_y/poly_size);

    //Check if there are more than 2 incident constrained edges.
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

    typename Gt::Construct_vector_2 vector = 
      Gt().construct_vector_2_object();

    // const Point_2& p = v->point();

    switch(num_constraints){
      // There aren't incident constraints.
      case 0: 
        if(cdt.is_inside_triangulation_cell(v,new_point)){
          typename CDT::Locate_type loc;
          int li;
          Face_handle fh = cdt.locate(new_point, loc, li);

          if( loc == CDT::VERTEX ){
            // Delete vertex on same position, or just return NULL?
            return Vector_2(v->point(),new_point);  
          }else{
            // Normal case
              return Vector_2(v->point(),new_point);
          }
        }
        else
          return CGAL::NULL_VECTOR;

      // This should be working as a Input Constraint
      case 1: return CGAL::NULL_VECTOR;

      // We need to check if the points is moving along the constraint.
      case 2: 
        if( CGAL::parallel(cdt.segment(parallel_constraint),new_segment)){
          
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

};

} // end namespace Mesh_2

} //namespace CGAL

#endif // CGAL_MESH_2_ODT_MOVE_H
