#ifndef CGAL_QT_CONSTRAINED_VORONOI_GRAPHICS_ITEM_H
#define CGAL_QT_CONSTRAINED_VORONOI_GRAPHICS_ITEM_H



#include <CGAL/Qt/GraphicsItem.h>
#include <CGAL/Qt/PainterOstream.h>
#include <CGAL/Qt/utility.h>

#include <QGraphicsScene>
#include <QGraphicsSceneMouseEvent>
#include <QPainter>
#include <QStyleOption>

#include <CGAL/intersection_2.h>
#include <CGAL/Polygon_2.h>

class QGraphicsSceneMouseEvent;


namespace CGAL {
namespace Qt {

template <typename DT>
class ConstrainedVoronoiGraphicsItem : public GraphicsItem
{
  typedef CGAL::Polygon_2<typename DT::Geom_traits, std::vector<typename DT::Point> >    Polygon_2;
public:
  ConstrainedVoronoiGraphicsItem(DT  * dt_);


  QRectF 
  boundingRect() const;
  
  void 
  paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *widget);
  
  void 
  modelChanged();

  const QPen& edgesPen() const
  {
    return edges_pen;
  }

  void setEdgesPen(const QPen& pen)
  {
    edges_pen = pen;
  }

private:
  DT * dt;
  QPen edges_pen;
};



template <typename DT>
ConstrainedVoronoiGraphicsItem<DT>::ConstrainedVoronoiGraphicsItem(DT * dt_)
  :  dt(dt_)
{
  setZValue(3);
}

template <typename DT>
QRectF 
ConstrainedVoronoiGraphicsItem<DT>::boundingRect() const
{
  QRectF rect = CGAL::Qt::viewportsBbox(scene());
  return rect;
}


template <typename DT>
void 
ConstrainedVoronoiGraphicsItem<DT>::paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget * /*w*/)
{
  QRectF rect = option->exposedRect;
  PainterOstream<typename DT::Geom_traits> pos(painter, rect);
  
  painter->setPen(edgesPen());
  for(typename DT::Finite_vertices_iterator vit = dt->finite_vertices_begin();
      vit != dt->finite_vertices_end();
      vit++){
    Polygon_2 poly = this->dt->dual(vit);
    for(unsigned int i=0; i<poly.size() ; i++){
      typename DT::Segment s = typename DT::Segment(poly[i],poly[(i+1)%poly.size()]);
      pos<<s;
    }
  }
}


template <typename T>
void 
ConstrainedVoronoiGraphicsItem<T>::modelChanged()
{
  update();
}

} // namespace Qt
} // namespace CGAL

#endif // CGAL_QT_CONSTRAINED_VORONOI_GRAPHICS_ITEM_H
