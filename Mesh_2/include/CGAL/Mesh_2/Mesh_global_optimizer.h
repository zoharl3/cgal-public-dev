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
  
  
template <typename CDT,
          typename MoveFunction>
class Mesh_global_optimizer
{  
  // Types
  typedef CDT                                              Tr;
  typedef typename Tr::Geom_traits                         Gt;
  
  typedef typename Tr::Point                               Point_2;
//  typedef typename Tr::Face_handle      Face_handle;
  typedef typename Tr::Vertex_handle                       Vertex_handle;
//  typedef typename Tr::Edge             Edge;
//  typedef typename Tr::Vertex           Vertex;
  
//  typedef typename Gt::FT                                  FT;
  typedef typename Gt::Vector_2                            Vector_2;
  
  typedef typename std::set<Vertex_handle>                 Vertex_set;
  typedef typename std::pair<Vertex_handle,Point_2>        Move;

  typedef std::vector<Move>                                Moves_vector;
  
//  typedef typename MoveFunction::Sizing_field Sizing_field;
      
public:
  /**
   * Constructor
   */
  Mesh_global_optimizer(CDT& cdt,
                        const MoveFunction);
  
  /**
   * Launch optimization process
   *
   * @param nb_interations maximum number of iterations
   */
  void operator()(int nb_iterations);

  /// Time accessors
  //void set_time_limit(double time) { time_limit_ = time; }
  //double time_limit() const { return time_limit_; }
  
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
   * update big_moves_ vector with new_sq_move value
   */
//  void update_big_moves(const FT& new_sq_move);

  /**
   * Updates mesh using moves of \c moves vector. Updates moving_vertices with
   * the new set of moving vertices after the move.
   */
  void update_mesh(const Moves_vector& moves,
                   Vertex_set& moving_vertices);

  /*bool is_time_limit_reached() const
  {
    return ( (this->time_limit() > 0) && (running_time_.time() > this->time_limit()) );
  }*/


private:
  // -----------------------------------
  // Private data
  // -----------------------------------
  CDT& cdt_;
//  FT sq_freeze_ratio_;
//  FT convergence_ratio_;
  MoveFunction move_function_;
//  Sizing_field sizing_field_;
//  double time_limit_;
//  CGAL::Timer running_time_;
  
//  typedef std::list<FT> FT_list;
//  FT_list big_moves_;
  
/*#ifdef CGAL_MESH_2_OPTIMIZER_VERBOSE
  mutable FT sum_moves_;
#endif*/
  
};
template <typename CDT, typename MoveFunction>
Mesh_global_optimizer<CDT, MoveFunction>::
Mesh_global_optimizer(CDT& cdt,
                        const MoveFunction move_function = MoveFunction())
: cdt_(cdt)
, move_function_(move_function)
{
}

template <typename CDT, typename MoveFunction>
void
Mesh_global_optimizer<CDT, MoveFunction>::
operator()(int nb_iterations)
{
  Vertex_set moving_vertices;

  for ( typename CDT::Finite_vertices_iterator vit = cdt_.finite_vertices_begin();
        vit != cdt_.finite_vertices_end();
        ++vit)
  {
    moving_vertices.insert(vit);
  }
  //collect_all_vertices(moving_vertices);

//  std::size_t initial_vertices_nb = moving_vertices.size();

  // Initialize big moves (stores the largest moves)
/*  big_moves_.clear();
  std::size_t big_moves_size_ = 
    (std::max)(std::size_t(1), std::size_t(moving_vertices.size()/500));
  big_moves_.resize(big_moves_size_,FT(0));
*/
  // Iterate
  int i = -1;
  while ( ++i < nb_iterations)// && ! is_time_limit_reached() )
  {
    // Compute move for each vertex
    Moves_vector moves = compute_moves(moving_vertices);
    //std::cout<<"MOVES SIZE: "<<moves.size()<<std::endl;

    // Stop if time_limit is reached
//    if ( is_time_limit_reached() )
//      break;
    
    // Update mesh with those moves
     // do : move_point(vertex, newpoint) for each pair    
    update_mesh(moves, moving_vertices); 
  }  
//  running_time_.stop();
  
}

template <typename CDT, typename MoveFunction>
typename Mesh_global_optimizer<CDT, MoveFunction>::Moves_vector
Mesh_global_optimizer<CDT, MoveFunction>::
compute_moves(const Vertex_set& moving_vertices)
{
  typename Gt::Construct_translated_point_2 translate =
    Gt().construct_translated_point_2_object();
  
  // Store new location of points which have to move
  Moves_vector moves;
  moves.reserve(moving_vertices.size());
  
  // reset worst_move list
//  std::fill(big_moves_.begin(),big_moves_.end(),FT(0));
  
  // Get move for each moving vertex
  int i=0;
  for ( typename Vertex_set::const_iterator vit = moving_vertices.begin() ;
       vit != moving_vertices.end() ;
       ++vit )
  {
    //std::cout<<++i<<": \n";
    Vector_2 move = compute_move(*vit);
    if ( CGAL::NULL_VECTOR != move )
    {
      //std::cout<<"BUENAS NOTICIAS"<<std::endl;
      Point_2 new_position = translate((*vit)->point(),move);
      moves.push_back(std::make_pair(*vit,new_position));
    }
    
    // Stop if time_limit_ is reached
    /*if ( is_time_limit_reached() )
      break;*/
  }   
  
  return moves;
}

template <typename CDT, typename MoveFunction>
typename Mesh_global_optimizer<CDT, MoveFunction>::Vector_2
Mesh_global_optimizer<CDT, MoveFunction>::
compute_move(const Vertex_handle& v)
{    
  /*typename Gt::Compute_squared_length_2 sq_length =
    Gt().compute_squared_length_2_object();*/
  
  //typename Gt::Construct_vector_2 vector =
    //Gt().construct_vector_2_object();
  
  /*typename Gt::Construct_translated_point_2 translate =
    Gt().construct_translated_point_2_object();*/
  
  // Get move from move function
  //Vector_2 move = move_function_(v, c2t2_, sizing_field_);
  Vector_2 move = move_function_(v, this->cdt_);
  
  // Project surface vertex
  /*if ( c2t2_.in_dimension(v) == 2 )
  {
    Point_2 new_position = translate(v->point(),move);

    //TODO: NEED TO DEFINE WHAT THIS MOVE DOES and how to get the vector_2.
    //move = vector(v->point(), helper_.project_on_surface(new_position,v));
  }*/
  
/*  FT local_sq_size = min_circumradius_sq_length(v);
  if ( FT(0) == local_sq_size )
    return CGAL::NULL_VECTOR;
  */
  //FT local_move_sq_length = sq_length(move) / local_sq_size;
  
  // Move point only if displacement is big enough w.r.t local size
  /*if ( local_move_sq_length < sq_freeze_ratio_ )
  {
    return CGAL::NULL_VECTOR;
  }*/
  
  // Update big moves
  //update_big_moves(local_move_sq_length);
  
  return move;
}
/*
template <typename CDT, typename MoveFunction>
void
Mesh_global_optimizer<CDT, MoveFunction>::
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
}*/

template <typename CDT, typename MoveFunction>
void
Mesh_global_optimizer<CDT, MoveFunction>::
update_mesh(const Moves_vector& moves,
            Vertex_set& moving_vertices)
{
  //std::set<Face_handle> outdated_faces;
  std::cout<<"moving point to its new location"<<std::endl;
  for( typename Moves_vector::const_iterator it = moves.begin() ;
        it != moves.end() ;
        ++it)
  {
    const Vertex_handle& v = it->first;
    const Point_2& new_position = it->second;

    // How to treat the sizing field?
    //move_point(v,new_position,outdated_faces);
    //std::cout<<"moving point to its new location: "<<new_position<<std::endl;
    cdt_.move(v,new_position);
  }
  moving_vertices.clear();
  // WHAT TO REALLY DO WITH MOVING_VERTICES?
    // Rebuild Delaunay?
} 

  
} // end namespace Mesh_2

} //namespace CGAL

#endif // CGAL_MESH_2_MESH_GLOBAL_OPTIMIZER_H
