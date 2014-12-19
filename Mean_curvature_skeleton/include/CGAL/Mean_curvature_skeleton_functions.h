// Copyright (c) 2014  GeometryFactory (France). All rights reserved.
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
// Author(s)     : Sebastien Loriot
//

#ifndef CGAL_MEAN_CURVATURE_SKELETON_FUNCTIONS_H
#define CGAL_MEAN_CURVATURE_SKELETON_FUNCTIONS_H

#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_items_with_id_3.h>
#include <CGAL/FaceGraph_to_Polyhedron_3.h>
#ifdef CGAL_EIGEN3_ENABLED
#include <CGAL/Eigen_solver_traits.h>  // for sparse linear system solver
#endif

#include <CGAL/Mean_curvature_skeleton.h>

namespace CGAL{

#if defined(DOXYGEN_RUNNING) || defined(CGAL_EIGEN3_ENABLED)
/// \ingroup PkgMeanCurvatureSkeleton3
/// @brief Extracts a medially centered curve skeleton for the mesh `hg`.
/// This function uses the class CGAL::Mean_curvature_flow_skeletonization with the default parameters.
/// This function is available if \ref thirdpartyEigen "Eigen" 3.2 (or greater) is available and `CGAL_EIGEN3_ENABLED` is defined.
/// @pre `hg` is a watertight triangulated surface mesh and has exactly one connected component.
/// @pre The specialization `boost::property_map<HalfedgeGraph, boost::vertex_point_t>::%type` and `get(vertex_point, hg)` are defined.
/// @pre The value type of `boost::property_map<HalfedgeGraph, boost::vertex_point_t>::%type` is a point type from a \cgal Kernel.
///
/// @tparam HalfedgeGraph
///         a model of `HalfedgeGraph`
/// @tparam Graph
///         an instantiation of `boost::adjacency_list` data structure for skeleton curve
/// @tparam VerticesPMap
///         a model of `ReadWritePropertyMap`</a>
///         with Graph::vertex_descriptor as key and
///         `std::vector< boost::graph_traits<HalfedgeGraph>::%vertex_descriptor >` as value type
/// @tparam PointPMap
///         a model of `ReadWritePropertyMap`</a>
///         with Graph::vertex_descriptor as key and
///         the value type of `boost::property_map<HalfedgeGraph, boost::vertex_point_t>::%type`
///         as value type.
///
/// @param hg
///        input mesh
/// @param skeleton
///        graph that will contain the skeleton of `hg`
/// @param skeleton_points
///        property map containing the location of the vertices of the graph `skeleton`
/// @param skeleton_to_hg_vertices property map associating a vertex `v` of the graph `skeleton`
///        to the set of vertices of `hg` corresponding to `v`.
/// \todo add an overload when the HalfedgeGraph is of the same type as the copy
/// \todo I need to tweak VerticesPMap to match the documentation
template <class HalfedgeGraph,
          class Graph,
          class PointPMap,
          class VerticesPMap>
void extract_mean_curvature_flow_skeleton(const HalfedgeGraph& hg,
                                          Graph& skeleton,
                                          PointPMap& skeleton_points,
                                          VerticesPMap& skeleton_to_hg_vertices)
{
  typedef typename boost::property_map<HalfedgeGraph, boost::vertex_point_t>::type HGPointPMap;
  typedef typename boost::property_traits<HGPointPMap>::value_type Point;
  typedef typename CGAL::Kernel_traits< Point >::Kernel K;
  typedef CGAL::Polyhedron_3<K,CGAL::Polyhedron_items_with_id_3> Polyhedron;

  // copy the input FaceGraph into a Polyhedron
  CGAL::FaceGraph_to_Polyhedron_3<HalfedgeGraph,
                                  HGPointPMap,
                                  typename Polyhedron::HalfedgeDS,
                                  false> modifier(hg, get(vertex_point, const_cast<HalfedgeGraph&>(hg)) );
  Polyhedron P;
  P.delegate(modifier);
  //init indices
  typedef typename boost::graph_traits<Polyhedron>::vertex_descriptor vertex_descriptor;
  typedef typename boost::graph_traits<Polyhedron>::halfedge_descriptor halfedge_descriptor;
  std::size_t i=0;
  BOOST_FOREACH( vertex_descriptor vd, vertices(P) )
    vd->id()=i++;
  i=0;
  BOOST_FOREACH( halfedge_descriptor hd, halfedges(P) )
    hd->id()=i++;

  typedef typename boost::property_map<Polyhedron, boost::vertex_index_t>::type VertexIndexMap;
  typedef typename boost::property_map<Polyhedron, boost::halfedge_index_t>::type HedgeIndexMap;
  typedef typename boost::property_map<Polyhedron, boost::vertex_point_t>::type PolyVertexPointMap;

  typedef CGAL::Eigen_solver_traits<
          Eigen::SparseLU<
          CGAL::Eigen_sparse_matrix<double>::EigenType,
          Eigen::COLAMDOrdering<int> > > SparseLinearAlgebraTraits_d;

  // extract the skeleton
  typedef CGAL::Mean_curvature_flow_skeletonization<Polyhedron,
                                                    VertexIndexMap,
                                                    HedgeIndexMap,
                                                    PolyVertexPointMap,
                                                    SparseLinearAlgebraTraits_d> Mcfskel;
  Mcfskel mcfs(P, get(boost::vertex_index,P), get(halfedge_index,P));

  mcfs.contract_until_convergence();
  mcfs.convert_to_skeleton(skeleton, skeleton_points, skeleton_to_hg_vertices);
}
#endif

}// end of namespace CGAL

#endif //CGAL_MEAN_CURVATURE_SKELETON_FUNCTIONS_H
