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
// Author(s)     : Stephane Tayeb
//
//******************************************************************************
// File Description : 
//******************************************************************************

#ifndef CGAL_MESH_2_MESH_GLOBAL_OPTIMIZER_H
#define CGAL_MESH_2_MESH_GLOBAL_OPTIMIZER_H

#ifdef CGAL_MESH_2_VERBOSE
  #define CGAL_MESH_2_OPTIMIZER_VERBOSE 
#endif

#include <CGAL/Timer.h>
#include <CGAL/Origin.h>

#include <vector>
#include <list>

#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>

namespace CGAL {

namespace Mesh_2 {
  
  
template <typename C2T2,
          typename MeshDomain,
          typename MoveFunction,
          typename Visitor_ = Null_global_optimizer_visitor<C3T3> >
class Mesh_global_optimizer
{  
  // Types
  typedef typename C2T2::Triangulation  Tr;
  typedef typename Tr::Geom_traits      Gt;
  
  typedef typename Tr::Point            Point_2;
  typedef typename Tr::Face_handle      Face_handle;
  typedef typename Tr::Vertex_handle    Vertex_handle;
  typedef typename Tr::Edge             Edge;
  typedef typename Tr::Vertex           Vertex;
  
  typedef typename Gt::FT               FT;
  typedef typename Gt::Vector_2         Vector_2;
  
  typedef typename std::vector<Face_handle>                 Face_vector;
  typedef typename std::vector<Vertex_handle>               Vertex_vector;
  typedef typename std::set<Vertex_handle>                  Vertex_set;
  typedef typename std::pair<Vertex_handle,Point_2>         Move;

  typedef std::vector<Move >   Moves_vector;
  
  typedef typename MoveFunction::Sizing_field Sizing_field;
  

  
  // Visitor class
  // Should define:
  //  - after_compute_moves()
  //  - after_move_points()
  //  - after_rebuild_restricted_delaunay()
  //  - end_of_iteration(int iteration_number)
  typedef Visitor_ Visitor;
    
public:
  /**
   * Constructor
   */
  Mesh_global_optimizer(C2T2& c2t2,
                        const MeshDomain& domain,
                        const FT& freeze_ratio,
                        const FT& convergence_ratio,
                        const MoveFunction move_function = MoveFunction());
  
  /// Time accessors
  void set_time_limit(double time) { time_limit_ = time; }
  double time_limit() const { return time_limit_; }
  

private:
  /**
   * Returns moves for vertices of set \c moving_vertices
   */
  Moves_vector compute_moves(const Vertex_set& moving_vertices);
  /**
   * Returns the move for vertex \c v
   */
  Vector_2 compute_move(const Vertex_handle& v);
  /**
   * Returns the minimum cicumradius length of faces incident to \c v
   */
  FT min_circumradius_sq_length(const Vertex_handle& v) const;
  /**
   * update big_moves_ vector with new_sq_move value
   */
  void update_big_moves(const FT& new_sq_move);

  bool is_time_limit_reached() const
  {
    return ( (time_limit() > 0) && (running_time_.time() > time_limit()) );      
  }


private:
  // -----------------------------------
  // Private data
  // -----------------------------------
  C2T2& c2t2_;
  Tr& tr_;
  const MeshDomain& domain_;
  FT sq_freeze_ratio_;
  FT convergence_ratio_;
  MoveFunction move_function_;
  Sizing_field sizing_field_;
  double time_limit_;
  CGAL::Timer running_time_;
  
  typedef std::list<FT> FT_list;
  FT_list big_moves_;
  
#ifdef CGAL_MESH_2_OPTIMIZER_VERBOSE
  mutable FT sum_moves_;
#endif



};

template <typename C2T2, typename Md, typename Mf, typename V_>
typename Mesh_global_optimizer<C2T2,Md,Mf,V_>::Moves_vector
Mesh_global_optimizer<C2T2,Md,Mf,V_>::
compute_moves(const Vertex_set& moving_vertices)
{
  typename Gt::Construct_translated_point_2 translate =
    Gt().construct_translated_point_2_object();
  
  // Store new location of points which have to move
  Moves_vector moves;
  moves.reserve(moving_vertices.size());
  
  // reset worst_move list
  std::fill(big_moves_.begin(),big_moves_.end(),FT(0));
  
  // Get move for each moving vertex
  for ( typename Vertex_set::const_iterator vit = moving_vertices.begin() ;
       vit != moving_vertices.end() ;
       ++vit )
  {
    Vector_2 move = compute_move(*vit);
    if ( CGAL::NULL_VECTOR != move )
    {
      Point_2 new_position = translate((*vit)->point(),move);
      moves.push_back(std::make_pair(*vit,new_position));
    }
    
    // Stop if time_limit_ is reached
    if ( is_time_limit_reached() )
      break;
  }   
  
  return moves;
}

template <typename C2T2, typename Md, typename Mf, typename V_>
typename Mesh_global_optimizer<C2T2,Md,Mf,V_>::Vector_2
Mesh_global_optimizer<C2T2,Md,Mf,V_>::
compute_move(const Vertex_handle& v)
{    
  typename Gt::Compute_squared_length_2 sq_length =
    Gt().compute_squared_length_2_object();
  
  typename Gt::Construct_vector_2 vector =
    Gt().construct_vector_2_object();
  
  typename Gt::Construct_translated_point_2 translate =
    Gt().construct_translated_point_2_object();
  
  // Get move from move function
  Vector_2 move = move_function_(v, c2t2_, sizing_field_);
  
  // Project surface vertex
  if ( c2t2_.in_dimension(v) == 2 )
  {
    Point_2 new_position = translate(v->point(),move);

    //TODO: NEED TO DEFINE WHAT THIS MOVE DOES and how to get the vector_2.
    //move = vector(v->point(), helper_.project_on_surface(new_position,v));
  }
  
  FT local_sq_size = min_circumradius_sq_length(v);
  if ( FT(0) == local_sq_size )
    return CGAL::NULL_VECTOR;
  
  FT local_move_sq_length = sq_length(move) / local_sq_size;
  
  // Move point only if displacement is big enough w.r.t local size
  if ( local_move_sq_length < sq_freeze_ratio_ )
  {
    return CGAL::NULL_VECTOR;
  }
  
  // Update big moves
  update_big_moves(local_move_sq_length);
  
  return move;
}

template <typename C2T2, typename Md, typename Mf, typename V_>
typename Mesh_global_optimizer<C2T2,Md,Mf,V_>::FT
Mesh_global_optimizer<C2T2,Md,Mf,V_>::
min_circumradius_sq_length(const Vertex_handle& v) const
{
  Face_vector incident_faces;
  incident_faces.reserve(64);
  tr_.incident_faces(v, std::back_inserter(incident_faces));
  
  // Get first face sq_circumradius_length 
  typename Face_vector::iterator cit = incident_faces.begin();
  //while ( incident_faces.end() != cit && !c2t2_.is_in_complex(*cit) ) { ++cit; }
  while ( incident_faces.end() != cits ) { ++cit; }
  
  // if vertex is isolated ...
  if ( incident_faces.end() == cit )
    return FT(0);
  
  // Initialize min
  FT min_sq_len = sq_circumradius_length(*cit++,v);
  
  // Find the minimum value
  for ( ; cit != incident_faces.end() ; ++cit )
  {
    //if ( !c2t2_.is_in_complex(*cit) )
      //continue;
    
    min_sq_len = (std::min)(min_sq_len,sq_circumradius_length(*cit,v));
  }
  
  return min_sq_len;
}

template <typename C2T2, typename Md, typename Mf, typename V_>
void
Mesh_global_optimizer<C2T2,Md,Mf,V_>::
update_big_moves(const FT& new_sq_move)
{  
  namespace bl = boost::lambda;
  
  if ( new_sq_move > big_moves_.back() )
  {
    // Remove last value
    big_moves_.pop_back();
    
    // Insert value at the right place
    typename FT_list::iterator pos = 
      std::find_if(big_moves_.begin(), big_moves_.end(), bl::_1 < new_sq_move );
    
    big_moves_.insert(pos, new_sq_move);
  }
}

  
} // end namespace Mesh_2

} //namespace CGAL

#endif // CGAL_MESH_2_MESH_GLOBAL_OPTIMIZER_H
