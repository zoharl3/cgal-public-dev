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
// File Description : Lloyd move function
//******************************************************************************

#ifndef CGAL_MESH_2_LLOYD_MOVE_H
#define CGAL_MESH_2_LLOYD_MOVE_H

#include <string>
#include <CGAL/centroid.h>
#include <CGAL/Polygon_2_algorithms.h>
#include <CGAL/intersections.h>
#include <CGAL/Mesh_2/Uniform_sizing_field.h>

namespace CGAL {

namespace Mesh_2 {

template <typename CDT,
  typename SizingField = Uniform_sizing_field<typename CDT::Triangulation> >
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

public:
  typedef SizingField Sizing_field;
  /**
   * @brief Return move to apply on \c v according to Lloyd optimization 
   * function
   */
  Vector_2 operator()(const Vertex_handle& v,
                       CDT& cdt,
                       const Sizing_field& sizing_field = Sizing_field() ) const
  {
    
    if (cdt.cell_is_infinite(v)){
      return CGAL::NULL_VECTOR;
    }
    if(v->input_constraint()){
      return CGAL::NULL_VECTOR;
    }

    Polygon_2 poly = cdt.dual(v);
    std::size_t poly_size = poly.size();

    if(poly_size==0){// This shouldn't happen
      return CGAL::NULL_VECTOR;
    }

    Point_2 new_point = CGAL::centroid(poly.vertices_begin(),
      poly.vertices_end(),
      CGAL::Dimension_tag<0>());

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
  static std::string name() { return std::string("Lloyd"); }
#endif
  
};

} // end namespace Mesh_2

} //namespace CGAL

#endif // CGAL_MESH_2_LLOYD_MOVE_H
