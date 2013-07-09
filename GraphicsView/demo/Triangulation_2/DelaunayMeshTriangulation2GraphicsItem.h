// Copyright (c) 2008  GeometryFactory Sarl (France).
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
// Author(s)     : Laurent Rineau <Laurent.Rineau@geometryfactory.com>

#ifndef CGAL_QT_DELAUNAY_MESH_TRIANGULATION2_GRAPHICS_ITEM_H
#define CGAL_QT_DELAUNAY_MESH_TRIANGULATION2_GRAPHICS_ITEM_H

#include <CGAL/Qt/ConstrainedTriangulationGraphicsItem.h>
#include <CGAL/Polygon_2.h>
#include <QBrush>

namespace CGAL {

namespace Qt {

template <typename T>
class DelaunayMeshTriangulation2GraphicsItem : public ConstrainedTriangulationGraphicsItem<T>
{
  typedef ConstrainedTriangulationGraphicsItem<T> Base;
  typedef CGAL::Polygon_2<typename T::Geom_traits, std::vector<typename T::Point> >    Polygon_2;

public:
  DelaunayMeshTriangulation2GraphicsItem(T  * t_)
    : Base(t_),
      visible_in_domain(true),
      in_domain_brush(::Qt::blue)
  {
  }


  const QBrush& facesInDomainBrush() const
  {
    return in_domain_brush;
  }

  void setFacesInDomainBrush(const QBrush& brush)
  {
    in_domain_brush = brush;
  }

  bool visibleFacesInDomain() const
  {
    return visible_in_domain;
  }

  void setVisibleFacesInDomain(const bool b)
  {
    visible_in_domain = b;
    this->update();
  }

protected:
  void drawAll(QPainter *painter);

  bool visible_in_domain;
  QBrush in_domain_brush;
};

template <typename T>
void 
DelaunayMeshTriangulation2GraphicsItem<T>::drawAll(QPainter *painter)
{
    this->painterostream = PainterOstream<typename T::Geom_traits>(painter);
    painter->setBrush(facesInDomainBrush());
    painter->setPen(::Qt::NoPen);
    
	//this->painterostream << this->t->dual(this->t->finite_vertices_begin());
  Polygon_2 poly = this->t->dual(this->t->finite_vertices_begin());
        for(typename Polygon_2::iterator fvit=poly.vertices_begin();
        fvit != poly.vertices_end();
        fvit++)
          this->painterostream << (*fvit)<<" ";
        //this->painterostream <<"\n";
  Base::drawAll(painter);
}

} // namespace Qt
} // namespace CGAL

#endif // CGAL_Q_DELAUNAY_MESH_TRIANGULATION_GRAPHICS_ITEM_H
