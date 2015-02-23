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

#ifndef CGAL_MCFSKEL_FIXED_EDGE_MAP_H
#define CGAL_MCFSKEL_FIXED_EDGE_MAP_H

/// @cond CGAL_DOCUMENT_INTERNAL

/** 
 * @file Fixed_edge_map.h
 * @brief This file contains the class to record if an edge is fixed.
 *
 * It is needed when using simplification package to do the edge collapse.
 */

namespace CGAL {
namespace internal {

// Map used to mark edges as fixed
#include <CGAL/Unique_hash_map.h>

#include <boost/graph/graph_traits.hpp>

//
// BGL property map which indicates whether an edge is border OR is marked as non-removable
//
template <class TriangleMesh>
class Fixed_edge_map : public boost::put_get_helper<bool, Fixed_edge_map<TriangleMesh> >
{
public:

  typedef boost::readable_property_map_tag                                 category;
  typedef bool                                                             value_type;
  typedef bool                                                             reference;
  typedef typename boost::graph_traits<TriangleMesh>::edge_descriptor     key_type;
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor halfedge_descriptor;

  Fixed_edge_map(const TriangleMesh& hg_) : hg(hg_), mFixed(false) {}

  reference operator[](key_type const& e) const
  {
    halfedge_descriptor h=halfedge(e,hg);
    return is_border(h,hg) || is_border(opposite(h,hg),hg) || is_fixed(h);
  }

  void set_is_fixed (halfedge_descriptor const& h, bool is)
  {
    mFixed[h] = is;
    mFixed[opposite(h,hg)] = is;
  }

  bool is_fixed(halfedge_descriptor const& h) const
  {
    return mFixed.is_defined(h) ? mFixed[h] : false;
  }

private:
  const TriangleMesh& hg;
  CGAL::Unique_hash_map<halfedge_descriptor, bool> mFixed;
};

} //namespace internal
} //namespace CGAL

/// @endcond

#endif //CGAL_MCFSKEL_FIXED_EDGE_MAP_H
