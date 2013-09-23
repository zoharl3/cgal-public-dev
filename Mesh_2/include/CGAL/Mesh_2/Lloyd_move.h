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
#include <vector>
#include <algorithm>
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
  
  typedef typename CDT::Vertex_handle         Vertex_handle;
  typedef typename CDT::Polygon_2             Polygon_2;
  typedef typename CDT::Point                 Point_2;
  typedef typename CDT::Edge_circulator       Edge_circulator;
  typedef typename CDT::Face_circulator       Face_circulator;
  typedef typename CDT::Edge                  Edge;
  typedef typename CDT::Face_handle           Face_handle;
  typedef typename CDT::Triangle              Triangle;

  typedef typename Gt::FT                     FT;
  typedef typename Gt::Vector_2               Vector_2;
  typedef typename Gt::Segment_2              Segment_2;

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
    // The infinite vertex
    if (cdt.cell_is_infinite(v)){
      return CGAL::NULL_VECTOR;
    }

    // If the vertex is an input_constraint
    if(v->input_constraint()){
      return CGAL::NULL_VECTOR;
    }

    Polygon_2 poly = cdt.dual(v);
    std::size_t poly_size = poly.size();

    CGAL_assertion(poly_size!=0);

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

    switch(num_constraints){
      case 0:{ 
        Vector_2 new_pos = calculate_vector(v, poly, cdt, sizing_field);

        //Point_2 new_position = translate(v->point(),new_pos);
        //std::cout<<"the new position is: "<<new_position<<std::endl;
        
        // Normal case
        return new_pos;
        
      }
      case 1: return CGAL::NULL_VECTOR;
      case 2: {
        Vector_2 new_pos = calculate_vector(v, poly, cdt, sizing_field);
        typename Gt::Construct_translated_point_2 translate =
          Gt().construct_translated_point_2_object();

        Point_2 new_position = translate(v->point(),new_pos);
        Segment_2 new_segment = Segment_2(new_position,v->point());

        if( CGAL::parallel(cdt.segment(parallel_constraint),new_segment)){
          //Face_handle fh = cdt.locate(new_point);
          return new_pos;
        }
        else
          return CGAL::NULL_VECTOR;
      }
      default: return CGAL::NULL_VECTOR;
    }
  }
  
#ifdef CGAL_MESH_2_OPTIMIZER_VERBOSE
  static std::string name() { return std::string("Lloyd"); }
#endif

private:
  /**
   * returns density_2d
   */
  template <typename Sizing_field>
  FT density_2d(const Point_2& p,
                const Face_handle& face,
                const Sizing_field& sizing_field) const
  {
    FT s = sizing_field(p,face);
    CGAL_assertion(!is_zero(s));

    // 1 / s^(d+2)
    return ( 1/(s*s*s*s) );
  }



  Vector_2 calculate_vector(const Vertex_handle& v,
    Polygon_2 polygon,
    CDT& cdt,
    const Sizing_field& sizing_field) const
  {

    Point_2 p = v->point();

    // Clean Polygon
    std::vector<Point_2> polygon_vector;
    for(typename Polygon_2::Vertex_iterator vit = polygon.vertices_begin();
        vit!=polygon.vertices_end(); ++vit){
      if( std::find(polygon_vector.begin(),polygon_vector.end(),*vit)
        == polygon_vector.end() )
        polygon_vector.push_back(*vit);
    }

    typename Gt::Construct_vector_2 vector =
      Gt().construct_vector_2_object();

    typename Gt::Compute_area_2 area = 
      Gt().compute_area_2_object();

    typename Gt::Construct_centroid_2 centroid = 
      Gt().construct_centroid_2_object();

    typename Gt::Construct_triangle_2 triangle =
      Gt().construct_triangle_2_object();


    // Move data
    Vector_2 move = CGAL::NULL_VECTOR;
    FT sum_masses(0);

    typename std::vector<Point_2>::iterator pit = polygon_vector.begin();

    Point_2 a = *pit;
    Point_2 b = a;

    while(pit != polygon_vector.end()){
      pit++;
      Point_2 c;
      if(pit == polygon_vector.end()){
        c = a;
      }else{
        c = *pit;
      }

      Triangle tri = triangle(p,b,c);
      Point_2 tri_centroid = centroid(tri);

      // Get face on the particular point
      typename CDT::Locate_type loc;
      int li;
      Face_handle fh = cdt.locate(c, loc, li);

      // Compute mass
      FT density = density_2d(tri_centroid, fh, sizing_field);
      FT abs_area = CGAL::abs(area(tri));
      FT mass = abs_area * density;

      move = move + mass * vector(p,tri_centroid);
      sum_masses += mass;

      b = c;
    }
    
    CGAL_assertion(sum_masses != 0.0);
    return move / sum_masses;
  }
  
};

} // end namespace Mesh_2

} //namespace CGAL

#endif // CGAL_MESH_2_LLOYD_MOVE_H
