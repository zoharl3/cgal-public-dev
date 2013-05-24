#include <QtCore/qglobal.h>

#include "Messages_interface.h"
#include "Scene_polyhedron_item_decorator.h"
#include "Polyhedron_demo_plugin_interface.h"
#include "ui_Smoothing_fairing_widget.h"
#include "Polyhedron_type.h"

//#include <CGAL/Fill_hole.h>
#include <CGAL/Fill_hole_Polyhedron_3.h>

#include <QTime>
#include <QAction>
#include <QMainWindow>
#include <QApplication>
#include <QDockWidget>
#include <QEvent>
#include <QKeyEvent>
#include <QMouseEvent>

#include <vector>
#include <algorithm>
#include <queue>

#include <QGLViewer/qglviewer.h>
#include "opengl_tools.h"

class Q_DECL_EXPORT Scene_polyhedron_selectable_item : public Scene_polyhedron_item_decorator
{
  Q_OBJECT
  typedef Polyhedron::Vertex_handle Vertex_handle;

public:
  Scene_polyhedron_selectable_item(Scene_polyhedron_item* poly_item, Ui::SmoothingFairing* ui_widget) 
    : Scene_polyhedron_item_decorator(poly_item), ui_widget(ui_widget)
  { 
    connect(poly_item, SIGNAL(selected_vertex(void*)), this, SLOT(vertex_has_been_selected(void*)));
    poly_item->enable_facets_picking(true);

    QGLViewer* viewer = *QGLViewer::QGLViewerPool().begin();
    viewer->installEventFilter(this);
  }

  void draw_edges() const {
    poly_item->direct_draw_edges();
    if(rendering_mode == Wireframe) {
      draw_ROI();
    }
  }  
  void draw() const {
    poly_item->draw();
    draw_ROI();
  }
  void draw_ROI() const {
    if(!selected_vertices.empty() && ui_widget->Show_ROI_check_box->isChecked()) {
      CGAL::GL::Color color;
      CGAL::GL::Point_size point_size; point_size.set_point_size(5);
      color.set_rgb_color(0, 1.f, 0);

      ::glBegin(GL_POINTS);
      for(std::set<Vertex_handle>::iterator 
        it = selected_vertices.begin(),
        end = selected_vertices.end();
      it != end; ++it)
      {
        const Kernel::Point_3& p = (*it)->point();
        ::glVertex3d(p.x(), p.y(), p.z());
      }
      ::glEnd();
    }
  }

public slots:
  void changed() {
    // do not use decorator function, which calls changed on poly_item which cause deletion of AABB
  }
  void vertex_has_been_selected(void* void_ptr) 
  {
    // get vertex descriptor
    Get_vertex_handle get_vertex_handle;  
    get_vertex_handle.vertex_ptr = static_cast<Polyhedron::Vertex*>(void_ptr);
    poly_item->polyhedron()->delegate(get_vertex_handle);
    Vertex_handle clicked_vertex = get_vertex_handle.vh;
    // use clicked_vertex, do what you want 
    bool is_insert = ui_widget->Insertion_radio_button->isChecked();
    int k_ring = ui_widget->Brush_size_spin_box->value();
    std::map<Vertex_handle, int> selection = extract_k_ring(*poly_item->polyhedron(), clicked_vertex, k_ring);
    bool any_change = false;
    if(is_insert) {
      for(std::map<Vertex_handle, int>::iterator it = selection.begin(); it != selection.end(); ++it) {
        any_change |= selected_vertices.insert(it->first).second;
      }
    }else {
      for(std::map<Vertex_handle, int>::iterator it = selection.begin(); it != selection.end(); ++it) {
        any_change |= selected_vertices.erase(it->first) != 0;
      }
    }
    if(any_change) { emit itemChanged(); }
  }

protected:
  bool eventFilter(QObject* /*target*/, QEvent *event)
  {
    // This filter is both filtering events from 'viewer' and 'main window'
    // key events
    if(event->type() == QEvent::KeyPress || event->type() == QEvent::KeyRelease)  {
      QKeyEvent *keyEvent = static_cast<QKeyEvent*>(event);
      Qt::KeyboardModifiers modifiers = keyEvent->modifiers();

      state.shift_pressing = modifiers.testFlag(Qt::ShiftModifier);
    }
    // mouse events
    if(event->type() == QEvent::MouseButtonPress || event->type() == QEvent::MouseButtonRelease) {
      QMouseEvent* mouse_event = static_cast<QMouseEvent*>(event);
      if(mouse_event->button() == Qt::LeftButton) {
        state.left_button_pressing = event->type() == QEvent::MouseButtonPress;
      }   
    }

    if(!poly_item->visible()) { return false; } // if not visible just update event state but don't do any action

    // use mouse move event for paint-like selection
    if(event->type() == QEvent::MouseMove &&
      (state.shift_pressing && state.left_button_pressing) )    
    { // paint with mouse move event 
      QMouseEvent* mouse_event = static_cast<QMouseEvent*>(event);
      QGLViewer* viewer = *QGLViewer::QGLViewerPool().begin();
      qglviewer::Camera* camera = viewer->camera();

      bool found = false;
      const qglviewer::Vec& point = camera->pointUnderPixel(mouse_event->pos(), found);
      if(found)
      {
        const qglviewer::Vec& orig = camera->position();
        const qglviewer::Vec& dir = point - orig;
        poly_item->select(orig.x, orig.y, orig.z, dir.x, dir.y, dir.z);
      }
    }//end MouseMove
    return false;
  }
  std::map<Vertex_handle, int> extract_k_ring(const Polyhedron &P, Vertex_handle v, int k)
  {
    std::map<Vertex_handle, int>  D;
    std::queue<Vertex_handle>     Q;
    Q.push(v); D[v] = 0;

    int dist_v;
    while( !Q.empty() && (dist_v = D[Q.front()]) < k ) {
      v = Q.front();
      Q.pop();

      Polyhedron::Halfedge_around_vertex_circulator e = v->vertex_begin();
      do
      {
        Vertex_handle new_v = e->opposite()->vertex();
        if(D.insert(std::make_pair(new_v, dist_v + 1)).second) {
          Q.push(new_v);
        }
      } while(++e != v->vertex_begin());
    }
    return D;
  }

// structs
  struct Mouse_keyboard_state
  {
    bool shift_pressing;
    bool left_button_pressing;

    Mouse_keyboard_state() 
      : shift_pressing(false), left_button_pressing(false) { }
  };

  struct Get_vertex_handle : public CGAL::Modifier_base<Polyhedron::HDS>
  {
    Polyhedron::Vertex* vertex_ptr;
    Vertex_handle vh;
    void operator()(Polyhedron::HDS& hds) {
      vh = hds.vertex_handle(vertex_ptr);
    }
  };

// members
  Ui::SmoothingFairing* ui_widget;
  Mouse_keyboard_state state;
public:
  std::set<Vertex_handle> selected_vertices;
};

class Polyhedron_demo_smoothing_fairing_plugin :
  public QObject,
  public Polyhedron_demo_plugin_interface
{
  Q_OBJECT
  Q_INTERFACES(Polyhedron_demo_plugin_interface)

public:
  bool applicable() const { return qobject_cast<Scene_polyhedron_item*>(scene->item(scene->mainSelectionIndex())); }
  void print_message(QString message) { messages->information(message);}
  QList<QAction*> actions() const { return QList<QAction*>() << actionSmoothingFairing; }

  void init(QMainWindow* mainWindow, Scene_interface* scene_interface, Messages_interface* m);
  Scene_polyhedron_selectable_item* convert_to_selectable_polyhedron(Scene_interface::Item_id i, Scene_polyhedron_item* poly_item);
  Scene_polyhedron_item* convert_to_plain_polyhedron(Scene_interface::Item_id i, Scene_polyhedron_selectable_item* selectable_poly);

public slots:
  void smoothing_fairing_action();
  void dock_widget_visibility_changed(bool visible);
  void on_Fair_button_clicked();
  void on_Set_all_vertices_button_clicked();
  void on_Clear_ROI_button_clicked();
  void on_Show_ROI_check_box_stateChanged(int state);
private:
  typedef Scene_interface::Item_id Item_id;

  QMainWindow* mw;
  Scene_interface* scene;
  Messages_interface* messages;
  QAction* actionSmoothingFairing;

  QDockWidget* dock_widget;
  Ui::SmoothingFairing* ui_widget;

}; // end Polyhedron_demo_smoothing_fairing_plugin

void Polyhedron_demo_smoothing_fairing_plugin::init(QMainWindow* mw,
                                      Scene_interface* scene_interface,
                                      Messages_interface* m)
{
  this->mw = mw;
  scene = scene_interface;
  messages = m;
  actionSmoothingFairing = new QAction(tr("Smoothing and Fairing"), mw);
  connect(actionSmoothingFairing, SIGNAL(triggered()),
          this, SLOT(smoothing_fairing_action()));

  dock_widget = new QDockWidget("Smoothing and Fairing", mw);
  dock_widget->setVisible(false);
  ui_widget = new Ui::SmoothingFairing();

  ui_widget->setupUi(dock_widget);
  mw->addDockWidget(Qt::LeftDockWidgetArea, dock_widget);

  connect(dock_widget, SIGNAL(visibilityChanged(bool)), this, SLOT(dock_widget_visibility_changed(bool)) );
  connect(ui_widget->Fair_button,  SIGNAL(clicked()), this, SLOT(on_Fair_button_clicked()));  
  connect(ui_widget->Set_all_vertices_button,  SIGNAL(clicked()), this, SLOT(on_Set_all_vertices_button_clicked()));  
  connect(ui_widget->Clear_ROI_button,  SIGNAL(clicked()), this, SLOT(on_Clear_ROI_button_clicked()));  
  connect(ui_widget->Show_ROI_check_box, SIGNAL(stateChanged(int)), this, SLOT(on_Show_ROI_check_box_stateChanged(int)));
}

void Polyhedron_demo_smoothing_fairing_plugin::smoothing_fairing_action(){
  if(dock_widget != NULL) { 
    dock_widget->show(); 
  }
}
void Polyhedron_demo_smoothing_fairing_plugin::on_Set_all_vertices_button_clicked() {
  int item_id = scene->mainSelectionIndex();
  Scene_polyhedron_selectable_item* poly_item = qobject_cast<Scene_polyhedron_selectable_item*>(scene->item(item_id));
  if(!poly_item) {
    print_message("Error: there is no selected polyhedron item!");
    return;
  }
  Polyhedron::Vertex_iterator vb(poly_item->polyhedron()->vertices_begin()), ve(poly_item->polyhedron()->vertices_end());
  for( ;vb != ve; ++vb) {
    poly_item->selected_vertices.insert(vb);
  }
}
void Polyhedron_demo_smoothing_fairing_plugin::on_Clear_ROI_button_clicked() {
  int item_id = scene->mainSelectionIndex();
  Scene_polyhedron_selectable_item* poly_item = qobject_cast<Scene_polyhedron_selectable_item*>(scene->item(item_id));
  if(!poly_item) {
    print_message("Error: there is no selected polyhedron item!");
    return;
  }
  poly_item->selected_vertices.clear();
}
void Polyhedron_demo_smoothing_fairing_plugin::on_Show_ROI_check_box_stateChanged(int /*state*/)
{
  for(Scene_interface::Item_id i = 0, end = scene->numberOfEntries(); i < end; ++i)
  {
    Scene_polyhedron_selectable_item* poly_item = qobject_cast<Scene_polyhedron_selectable_item*>(scene->item(i));
    if(!poly_item) { continue; }

    scene->itemChanged(poly_item);  // just for redraw   
  }  
}

void Polyhedron_demo_smoothing_fairing_plugin::on_Fair_button_clicked() {

  int item_id = scene->mainSelectionIndex();
  Scene_polyhedron_selectable_item* poly_item = qobject_cast<Scene_polyhedron_selectable_item*>(scene->item(item_id));
  if(!poly_item) {
    print_message("Error: there is no selected polyhedron item!");
    return;
  }

  if(ui_widget->Scale_dependent_weight_radio_button->isChecked())
    CGAL::fair(*poly_item->polyhedron(), poly_item->selected_vertices, 
      CGAL::internal::Fairing_weight_selector<Polyhedron, CGAL::SCALE_DEPENDENT_WEIGHTING>::weight_calculator());
  if(ui_widget->Uniform_weight_radio_button->isChecked())
    CGAL::fair(*poly_item->polyhedron(), poly_item->selected_vertices,
      CGAL::internal::Fairing_weight_selector<Polyhedron, CGAL::UNIFORM_WEIGHTING>::weight_calculator());
  else
    CGAL::fair(*poly_item->polyhedron(), poly_item->selected_vertices,
      CGAL::internal::Fairing_weight_selector<Polyhedron, CGAL::COTANGENT_WEIGHTING>::weight_calculator());
  poly_item->changed();
}

Scene_polyhedron_selectable_item* 
  Polyhedron_demo_smoothing_fairing_plugin::convert_to_selectable_polyhedron(Item_id i,
  Scene_polyhedron_item* poly_item)
{
  QString poly_item_name = poly_item->name();
  Scene_polyhedron_selectable_item* selectable_poly = new Scene_polyhedron_selectable_item(poly_item, ui_widget);
  selectable_poly->setColor(poly_item->color());
  selectable_poly->setName(QString("%1 (selectable)").arg(poly_item->name()));

  mw->installEventFilter(selectable_poly); // filter mainwindows events for key(pressed/released)

  scene->replaceItem(i, selectable_poly);
  return selectable_poly;
}

Scene_polyhedron_item* 
  Polyhedron_demo_smoothing_fairing_plugin::convert_to_plain_polyhedron(Item_id i,
  Scene_polyhedron_selectable_item* selectable_poly) 
{
  Scene_polyhedron_item* poly_item = selectable_poly->to_polyhedron_item();
  scene->replaceItem(i, poly_item);
  delete selectable_poly;
  return poly_item;
}

void Polyhedron_demo_smoothing_fairing_plugin::dock_widget_visibility_changed(bool visible)
{
  for(Scene_interface::Item_id i = 0, end = scene->numberOfEntries();
    i < end; ++i)
  {
    Scene_polyhedron_item* poly_item = qobject_cast<Scene_polyhedron_item*>(scene->item(i));
    Scene_polyhedron_selectable_item* selectable_item = qobject_cast<Scene_polyhedron_selectable_item*>(scene->item(i));

    if(visible && poly_item) {
      convert_to_selectable_polyhedron(i, poly_item);
    } else if(!visible && selectable_item) {
      convert_to_plain_polyhedron(i, selectable_item);
    }
  }
}

Q_EXPORT_PLUGIN2(Polyhedron_demo_smoothing_fairing_plugin, Polyhedron_demo_smoothing_fairing_plugin)

#include "Polyhedron_demo_smoothing_fairing_plugin.moc"
