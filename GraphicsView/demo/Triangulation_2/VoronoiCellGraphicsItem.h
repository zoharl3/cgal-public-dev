#ifndef CGAL_QT_VORONOI_CELL_GRAPHICS_ITEM_H
#define CGAL_QT_VORONOI_CELL_GRAPHICS_ITEM_H



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
class VoronoiCellGraphicsItem : public GraphicsItem
{
  typedef CGAL::Polygon_2<typename DT::Geom_traits, std::vector<typename DT::Point> >    Polygon_2;
public:
  VoronoiCellGraphicsItem(DT  * dt_);


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
VoronoiCellGraphicsItem<DT>::VoronoiCellGraphicsItem(DT * dt_)
  :  dt(dt_)
{
  setZValue(3);
}

template <typename DT>
QRectF 
VoronoiCellGraphicsItem<DT>::boundingRect() const
{
  QRectF rect = CGAL::Qt::viewportsBbox(scene());
  return rect;
}


template <typename DT>
void 
VoronoiCellGraphicsItem<DT>::paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget * /*w*/)
{
  QRectF rect = option->exposedRect;
  PainterOstream<typename DT::Geom_traits> pos(painter, rect);
  
  painter->setPen(edgesPen());
  for(typename DT::Finite_edges_iterator eit = dt->finite_edges_begin();
      eit != dt->finite_edges_end();
      eit++){
    Polygon_2 poly = this->dt->dual(this->dt->finite_vertices_begin());
    //CGAL::Object o = dt->dual(eit);
    typename DT::Segment s;
    typename DT::Geom_traits::Ray_2 r;
    typename DT::Geom_traits::Line_2 l;
    //if(CGAL::assign(s,o)){
      //pos << s;
    //} else if(CGAL::assign(r,o)) {
//      pos << r;
    //}else if(CGAL::assign(l,o)) {
//      pos << l;
    //} 
  }
}


template <typename T>
void 
VoronoiCellGraphicsItem<T>::modelChanged()
{
  update();
}

} // namespace Qt
} // namespace CGAL

#endif // CGAL_QT_VORONOI_CELL_GRAPHICS_ITEM_H
