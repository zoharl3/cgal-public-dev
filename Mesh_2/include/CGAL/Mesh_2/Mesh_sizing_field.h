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
// Author(s)     : Raul Gallegos
//
//******************************************************************************
// File Description : Defines a sizing field stored into an external
// mesh triangulation
//******************************************************************************

#ifndef CGAL_MESH_2_MESH_SIZING_FIELD_H
#define CGAL_MESH_2_MESH_SIZING_FIELD_H

namespace CGAL {

namespace Mesh_2
{
  
/**
 * @class Mesh_sizing_field
 */
template <typename Tr, bool Need_vertex_update = true>
class Mesh_sizing_field
{
  // Types
  typedef typename Tr::Geom_traits        Gt;
  typedef typename Tr::Point              Point_2;
  typedef typename Gt::FT                 FT;

  typedef typename Tr::Vertex_handle      Vertex_handle;
  typedef typename Tr::Face_handle        Face_handle;
  
public:
  // update vertices of mesh triangulation ? 
  static const bool is_vertex_update_needed = Need_vertex_update;
  
public:
  /**
   * Constructor
   */
  Mesh_sizing_field(Tr& tr);
  
  /**
   * Fill sizing field, using size associated to point in \c value_map
   */
  void fill(const std::map<Point_2, FT>& value_map);
  
  /**
   * Returns size at point \c p.
   */
  FT operator()(const Point_2& p) const
  { return this->operator()(p,last_face_); }
  
  /**
   * Returns size at point \c p, using \c v to accelerate \c p location
   * in triangulation
   */
  FT operator()(const Point_2& p, const Vertex_handle& v) const
  //{ return this->operator(p,v->); }
  { return this->operator()(p, v->face());}
  
  /**
   * Returns size at point \c p.
   */
  FT operator()(const Point_2& p, const Face_handle& c) const;
  
  /**
   * Returns size at point \c p. Assumes that p is the centroid of c.
   */
  //FT operator()(const Point_2& p, const std::pair<Cell_handle,bool>& c) const;
  
private:
  /**
   * Returns size at point \c p, by interpolation into triangle.
   */
  FT interpolate_normal(const Point_2& p,
                                  const Face_handle& face) const;

  /**
   * Returns size at point \c p, by interpolation into triangle with the infinite vertex.
   */
  FT interpolate_infinite(const Point_2& p,
                                  const Face_handle& face) const;
  
  /**
   * Returns size at point \c p, by interpolation into facet (\c cell is assumed
   * to be an infinite cell).
   */
  //FT interpolate_on_facet_vertices(const Point_2& p,
    //                               const Cell_handle& cell) const;
  
private:
  /// The triangulation
  Tr& tr_;
  /// A cell that is used to accelerate location queries
  mutable Face_handle last_face_;
};
  
  
  
template <typename Tr, bool B>
Mesh_sizing_field<Tr,B>::
Mesh_sizing_field(Tr& tr)
  : tr_(tr)
  , last_face_()
{
}  
  
  
template <typename Tr, bool B>
void
Mesh_sizing_field<Tr,B>::
fill(const std::map<Point_2, FT>& value_map)
{
  typedef typename Tr::Finite_vertices_iterator  Fvi;
  
  for ( Fvi vit = tr_.finite_vertices_begin() ;
        vit != tr_.finite_vertices_end() ;
        ++ vit )
  {
    typename std::map<Point_2, FT>::const_iterator find_result = 
      value_map.find(vit->point());
    
    if ( find_result != value_map.end() )
    {
      vit->set_meshing_info(find_result->second);
    }
    else
    {
      CGAL_assertion(false);
      vit->set_meshing_info(FT(0));
    }
  }
}  

template <typename Tr, bool B>
typename Mesh_sizing_field<Tr,B>::FT
Mesh_sizing_field<Tr,B>::
operator()(const Point_2& p, const Face_handle& f) const  
{  
  const Face_handle face = tr_.locate(p,f);
  last_face_ = face;
  
  if( !tr_.is_infinite(face) )
    return interpolate_normal(p,face);
  else
    return interpolate_infinite(p,face);
}

/*
template <typename Tr, bool B>
typename Mesh_sizing_field<Tr,B>::FT
Mesh_sizing_field<Tr,B>::
operator()(const Point_2&, const std::pair<Cell_handle,bool>& c) const
{
  // Assumes that p is the centroid of c
  const Cell_handle& cell = c.first;
  
  // Interpolate value using tet vertices values
  const FT& va = cell->vertex(0)->meshing_info();
  const FT& vb = cell->vertex(1)->meshing_info();
  const FT& vc = cell->vertex(2)->meshing_info();
  const FT& vd = cell->vertex(3)->meshing_info();
  
  return ( (va+vb+vc+vd)/4 );
}*/

//TODO: check this method functionality.  
template <typename Tr, bool B>
typename Mesh_sizing_field<Tr,B>::FT
Mesh_sizing_field<Tr,B>::
interpolate_normal(const Point_2& p, const Face_handle& face) const
{
  typename Gt::Compute_area_2 area =
    Gt().compute_area_2_object();
  
  // Interpolate value using triangle vertices values
  const FT& va = face->vertex(0)->meshing_info();
  const FT& vb = face->vertex(1)->meshing_info();
  const FT& vc = face->vertex(2)->meshing_info();
  
  const Point_2& a = face->vertex(0)->point();
  const Point_2& b = face->vertex(1)->point();
  const Point_2& c = face->vertex(2)->point();
  
  const FT abp = CGAL::abs(area(a,b,p));
  const FT acp = CGAL::abs(area(a,c,p));
  const FT bcp = CGAL::abs(area(b,c,p));
  
  // If area is 0, then compute the average value
  if ( is_zero(abp+acp+bcp) )
    return (va+vb+vc)/3.;
  
  return ( (abp*vc + acp*vb + bcp*va) / (abp+acp+bcp) );
}

//TODO: check this method functionality.  
template <typename Tr, bool B>
typename Mesh_sizing_field<Tr,B>::FT
Mesh_sizing_field<Tr,B>::
interpolate_infinite(const Point_2& p, const Face_handle& face) const
{
  /*typename Gt::Compute_area_2 area =
    Gt().compute_area_2_object();

  // Find infinite vertex and put it in k0
  int k0 = 0;
  int k1 = 1;
  int k2 = 2;

  if ( tr_.is_infinite(face->vertex(1)) )
    std::swap(k0,k1);
  if ( tr_.is_infinite(face->vertex(2)) )
    std::swap(k0,k2);
  
  // Interpolate value using triangle vertices values
  const FT& va = face->vertex(k1)->meshing_info();
  const FT& vb = face->vertex(k2)->meshing_info();
  
  const Point_2& a = face->vertex(k1)->point();
  const Point_2& b = face->vertex(k2)->point();
  
  const FT abp = CGAL::abs(area(a,b,p));
  const FT acp = CGAL::abs(area(a,c,p));
  const FT bcp = CGAL::abs(area(b,c,p));
  
  // If area is 0, then compute the average value
  if ( is_zero(abp+acp+bcp) )
    return (va+vb+vc)/3.;
  
  return ( (abp*vc + acp*vb + bcp*va) / (abp+acp+bcp) );
  */
  std::cout<<"SHOULD THIS HAPPEN?"<<std::endl;
  FT ft(0);
  return ft;
}
  
  
/*  
template <typename Tr, bool B>
typename Mesh_sizing_field<Tr,B>::FT
Mesh_sizing_field<Tr,B>::
interpolate_on_facet_vertices(const Point_2& p, const Cell_handle& cell) const
{
  typename Gt::Compute_area_3 area =
    Gt().compute_area_3_object();
  
  // Find infinite vertex and put it in k0
  int k0 = 0;
  int k1 = 1;
  int k2 = 2;
  int k3 = 3;
  
  if ( tr_.is_infinite(cell->vertex(1)) )
    std::swap(k0,k1);
  if ( tr_.is_infinite(cell->vertex(2)) )
    std::swap(k0,k2);
  if ( tr_.is_infinite(cell->vertex(3)) )
    std::swap(k0,k3);
  
  // Interpolate value using tet vertices values
  const FT& va = cell->vertex(k1)->meshing_info();
  const FT& vb = cell->vertex(k2)->meshing_info();
  const FT& vc = cell->vertex(k3)->meshing_info();
  
  const Point_2& a = cell->vertex(k1)->point();
  const Point_2& b = cell->vertex(k2)->point();
  const Point_2& c = cell->vertex(k3)->point();
  
  const FT abp = area(a,b,p);
  const FT acp = area(a,c,p);
  const FT bcp = area(b,c,p);
  
  CGAL_assertion(abp >= 0);
  CGAL_assertion(acp >= 0);
  CGAL_assertion(bcp >= 0);
  
  // If area is 0, then compute the average value
  if ( is_zero(abp+acp+bcp) )
    return (va+vb+vc)/3.;
  
  return ( (abp*vc + acp*vb + bcp*va ) / (abp+acp+bcp) );
}*/
  
} // end namespace Mesh_2


} //namespace CGAL

#endif // CGAL_MESH_2_MESH_SIZING_FIELD_H
