// Copyright (c) 2013  GeometryFactory (France). All rights reserved.
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
// Author(s)     : Xiang Gao <gaox@ethz.ch>
//

#ifndef CGAL_MCFSKEL_UTILITY_H
#define CGAL_MCFSKEL_UTILITY_H

/// @cond CGAL_DOCUMENT_INTERNAL

/** 
 * @file Utility.h
 * @brief This file contains some helper functions like splitting an edge at a 
 * given point.
 */

#include <boost/graph/graph_traits.hpp>

namespace CGAL {
namespace internal {
  
/**
* Split the edge
* @param polyhedron the mesh containing the given edge
* @param ei the edge to be split
* @param pn the position of the new vertex created by the split
*/
template<class HalfedgeGraph, class HalfedgeGraphPointPMap>
typename boost::graph_traits<HalfedgeGraph>::halfedge_descriptor
mesh_split(HalfedgeGraph& polyhedron, 
           HalfedgeGraphPointPMap& hg_point_pmap,
           typename boost::graph_traits<HalfedgeGraph>::halfedge_descriptor ei,
           typename HalfedgeGraph::Traits::Point_3 pn)
{
  typedef typename boost::graph_traits<HalfedgeGraph>::halfedge_descriptor            halfedge_descriptor;

  
  // halfedge_descriptor en = Euler::split_edge(ei, polyhedron); // there is an issue in this function for now use the polyhedron version in the meantime
  halfedge_descriptor en = polyhedron.split_edge(ei);
  boost::put(hg_point_pmap, en->vertex(), pn);
  Euler::split_face(en, ei->next(), polyhedron);

  en->id() = -1;
  en->opposite()->id() = -1;
  ei->id() = -1;
  ei->opposite()->id() = -1;
  en->next()->id() = -1;
  en->next()->opposite()->id() = -1;
  en->next()->next()->id() = -1;
  ei->next()->id() = -1;
  halfedge_descriptor ej = opposite(en, polyhedron);
  if (!(ej->is_border()))
  {
    Euler::split_face(ei->opposite(), ej->next(), polyhedron);
    ej->next()->id() = -1;
    halfedge_descriptor ei_op_next = ei->opposite()->next();
    ei_op_next->id() = -1;
    ei_op_next->opposite()->id() = -1;
    ei_op_next->next()->id() = -1;
  }

  return en;
}

template<class Vertex, class Kernel>
double get_triangle_area(typename Kernel::Point_3 p1,
                         typename Kernel::Point_3 p2,
                         typename Kernel::Point_3 p3)
{
  typedef typename Kernel::Vector_3 Vector;
  Vector v12(p1, p2);
  Vector v13(p1, p3);
  return sqrtf(cross_product(v12, v13).squared_length()) * 0.5;
}

template<class HalfedgeGraph, class HalfedgeGraphPointPMap>
double get_surface_area(HalfedgeGraph& polyhedron, HalfedgeGraphPointPMap& hg_point_pmap)
{
  typedef typename HalfedgeGraph::Traits                                  Kernel;
  typedef typename Kernel::Point_3                                        Point;
  typedef typename HalfedgeGraph::Facet_iterator                          Facet_iterator;
  typedef typename HalfedgeGraph::Halfedge_around_facet_circulator        Halfedge_facet_circulator;
  typedef typename boost::graph_traits<HalfedgeGraph>::vertex_descriptor	vertex_descriptor;

  double total_area = 0;
  for (Facet_iterator i = polyhedron.facets_begin(); i != polyhedron.facets_end(); ++i)
  {
    Halfedge_facet_circulator j = i->facet_begin();
    vertex_descriptor v1 = j->vertex();
    ++j;
    vertex_descriptor v2 = j->vertex();
    ++j;
    vertex_descriptor v3 = j->vertex();
    Point p1 = boost::get(hg_point_pmap, v1);
    Point p2 = boost::get(hg_point_pmap, v2);
    Point p3 = boost::get(hg_point_pmap, v3);
    total_area += internal::get_triangle_area<vertex_descriptor, Kernel>(p1, p2, p3);
  }
  return total_area;
}

} //namespace internal
} //namespace CGAL

/// @endcond

#endif //CGAL_MCFSKEL_UTILITY_H
