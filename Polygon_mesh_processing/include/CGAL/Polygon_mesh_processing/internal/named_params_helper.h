// Copyright (c) 2013 GeometryFactory (France).
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
// Author(s)     : Jane Tournois

#ifndef CGAL_NAMED_PARAMETERS_HELPERS_H
#define CGAL_NAMED_PARAMETERS_HELPERS_H

#include <CGAL/Kernel_traits.h>
#include <CGAL/Polygon_mesh_processing/internal/named_function_params.h>

#include <CGAL/property_map.h>

// shortcut for accessing the value type of the property map
template <class Graph, class Property>
class property_map_value {
  typedef typename boost::property_map<Graph, Property>::const_type PMap;
public:
  typedef typename boost::property_traits<PMap>::value_type type;
};

template<typename PolygonMesh>
class GetK
{
  typedef typename property_map_value<PolygonMesh,
    boost::vertex_point_t>::type
    Point;
public:
  typedef typename CGAL::Kernel_traits<Point>::Kernel Kernel;
};

template<typename PolygonMesh, typename NamedParameters>
class GetKernel
{
  typedef typename GetK<PolygonMesh>::Kernel DefaultKernel;
public:
  typedef typename boost::lookup_named_param_def <
    CGAL::geom_traits_t,
    NamedParameters,
    DefaultKernel
  > ::type  Kernel;
};

template<typename PolygonMesh, typename NamedParameters>
class GetVertexPointMap
{
  typedef typename boost::property_map<PolygonMesh, boost::vertex_point_t>::const_type DefaultVPMap;
public:
  typedef typename boost::lookup_named_param_def<
    boost::vertex_point_t,
    NamedParameters,
    DefaultVPMap
  > ::type  type;

};

#endif //CGAL_NAMED_PARAMETERS_HELPERS_H

