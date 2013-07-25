#include "Polyhedron_demo_plugin_helper.h"
#include "Polyhedron_demo_plugin_interface.h"
#include "ui_Mean_curvature_flow_skeleton_plugin.h"
#include "Scene_polyhedron_item.h"
#include "Scene_points_with_normal_item.h"
#include "Scene_polylines_item.h"
#include "Scene.h"

#include "Polyhedron_type.h"
#include "MainWindow.h"
#include "ui_MainWindow.h"

#include <QApplication>
#include <QMainWindow>
#include <QInputDialog>
#include <QTime>
#include <QMessageBox>

#include <Eigen/Sparse>

#include <boost/property_map/property_map.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>

#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_items_with_id_3.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Eigen_solver_traits.h>
#include <CGAL/Mean_curvature_skeleton.h>
#include <CGAL/iterator.h>
#include <CGAL/internal/corefinement/Polyhedron_subset_extraction.h>

#include <queue>

template<class PolyhedronWithId, class KeyType>
struct Polyhedron_with_id_property_map
    : public boost::put_get_helper<std::size_t&,
             Polyhedron_with_id_property_map<PolyhedronWithId, KeyType> >
{
public:
    typedef KeyType      key_type;
    typedef std::size_t  value_type;
    typedef value_type&  reference;
    typedef boost::lvalue_property_map_tag category;

    reference operator[](key_type key) const { return key->id(); }
};

typedef boost::graph_traits<Polyhedron>::vertex_descriptor          vertex_descriptor;
typedef boost::graph_traits<Polyhedron>::vertex_iterator            vertex_iterator;
typedef boost::graph_traits<Polyhedron>::edge_descriptor            edge_descriptor;
typedef Polyhedron::Facet_iterator                                  Facet_iterator;
typedef Polyhedron::Halfedge_around_facet_circulator                Halfedge_facet_circulator;

typedef Polyhedron_with_id_property_map<Polyhedron, vertex_descriptor> Vertex_index_map; // use id field of vertices
typedef Polyhedron_with_id_property_map<Polyhedron, edge_descriptor>   Edge_index_map;   // use id field of edges

typedef CGAL::Eigen_solver_traits<Eigen::SimplicialLDLT<CGAL::Eigen_sparse_matrix<double>::EigenType> > Sparse_linear_solver;

typedef CGAL::Mean_curvature_skeleton<Polyhedron, Sparse_linear_solver, Vertex_index_map, Edge_index_map> Mean_curvature_skeleton;

typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS> Graph;

typedef boost::graph_traits<Graph>::in_edge_iterator                in_edge_iter;
typedef boost::graph_traits<Graph>::out_edge_iterator               out_edge_iter;
typedef boost::graph_traits<Graph>::edge_iterator                   edge_iter;
typedef boost::graph_traits<Graph>::edge_descriptor                 edge_desc;

typedef Polyhedron::Traits         Kernel;
typedef Kernel::Point_3            Point;

class Polyhedron_demo_mean_curvature_flow_skeleton_plugin :
  public QObject,
  public Polyhedron_demo_plugin_helper
{
  Q_OBJECT
  Q_INTERFACES(Polyhedron_demo_plugin_interface)
  QAction* actionMCFSkeleton;
  QAction* actionConvert_to_skeleton;
  QAction* actionConvert_to_medial_skeleton;

public:
  // used by Polyhedron_demo_plugin_helper
  QStringList actionsNames() const {
    return QStringList() << "actionMCFSkeleton" << "actionConvert_to_skeleton"
                         << "actionConvert_to_medial_skeleton";
  }

  void init(QMainWindow* mainWindow, Scene_interface* scene_interface) {
    mcs = NULL;
    dockWidget = NULL;
    ui = NULL;

    actionMCFSkeleton = new QAction(tr("Mean Curvature Skeleton"), mainWindow);
    actionMCFSkeleton->setObjectName("actionMCFSkeleton");

    actionConvert_to_skeleton = new QAction(tr("Extract Skeleton"), mainWindow);
    actionConvert_to_skeleton->setObjectName("actionConvert_to_skeleton");

    actionConvert_to_medial_skeleton = new QAction(tr("Extract Medial Skeleton"), mainWindow);
    actionConvert_to_medial_skeleton->setObjectName("actionConvert_to_medial_skeleton");

    Polyhedron_demo_plugin_helper::init(mainWindow, scene_interface);
  }

  QList<QAction*> actions() const {
    return QList<QAction*>() << actionMCFSkeleton << actionConvert_to_skeleton
                             << actionConvert_to_medial_skeleton;
  }

  bool applicable() const {
    return qobject_cast<Scene_polyhedron_item*>(scene->item(scene->mainSelectionIndex()));
  }

  void init_ui(double diag) {
    ui->omega_L->setValue(1);
    ui->omega_L->setSingleStep(0.1);
    ui->omega_H->setValue(0.1);
    ui->omega_H->setSingleStep(0.1);
    ui->omega_P->setValue(0.2);
    ui->omega_P->setSingleStep(0.1);
    ui->edgelength_TH->setDecimals(7);
    ui->edgelength_TH->setValue(0.002 * diag);
    ui->edgelength_TH->setSingleStep(0.0000001);
    ui->volume_TH->setDecimals(6);
    ui->volume_TH->setValue(1e-04);
    ui->volume_TH->setSingleStep(0.000001);
    ui->area_TH->setDecimals(7);
    ui->area_TH->setValue(1e-4);
    ui->area_TH->setSingleStep(1e-5);
    ui->is_medially_centered->setChecked(false);

    ui->label_omega_L->setToolTip(QString("omega_L / omega_H controls the velocity of movement and approximation quality"));
    ui->label_omega_H->setToolTip(QString("omega_L / omega_H controls the velocity of movement and approximation quality"));
    ui->label_omega_P->setToolTip(QString("omega_L / omega_P controls the smoothness of the medial approximation"));
    ui->label_volume_TH->setToolTip(QString("Run to converge will stop when (current volume / original volume) < volume_TH"));
    ui->pushButton_contract->setToolTip(QString("contract mesh based on mean curvature flow"));
    ui->pushButton_collapse->setToolTip(QString("collapse short edges"));
    ui->pushButton_split->setToolTip(QString("split obtuse triangles"));
    ui->pushButton_degeneracy->setToolTip(QString("fix degenerate points"));
    ui->pushButton_skeletonize->setToolTip(QString("Turn mesh to a skeleton curve"));
    ui->pushButton_run->setToolTip(QString("run one iteration of contract, collapse, split, detect degeneracy"));
    ui->pushButton_converge->setToolTip(QString("iteratively contract the mesh until convergence"));

    // only for debugging
    ui->pushButton_voronoi->setVisible(false);
  }

  bool check_item_index(int index) {
    if (index < 0)
    {
      QMessageBox msgBox;
      msgBox.setText("Please select an item first :)");
      msgBox.exec();
      return false;
    }
    return true;
  }

  bool is_mesh_valid(Polyhedron *pMesh) {
    if (!pMesh->is_closed())
    {
      QMessageBox msgBox;
      msgBox.setText("The mesh is not closed.");
      msgBox.exec();
      return false;
    }
    if (!pMesh->is_pure_triangle())
    {
      QMessageBox msgBox;
      msgBox.setText("The mesh is not a pure triangle mesh.");
      msgBox.exec();
      return false;
    }

    std::size_t num_component;
    CGAL::Counting_output_iterator output_it(&num_component);
    CGAL::internal::extract_connected_components(*pMesh, output_it);
    ++output_it;
    if (num_component != 1)
    {
      QMessageBox msgBox;
      QString str = QString("The mesh is not a single closed mesh.\n It has %1 components.").arg(num_component);
      msgBox.setText(str);
      msgBox.exec();
      return false;
    }
    return true;
  }

  // check if the Mean_curvature_skeleton exists
  // or has the same polyheron item
  // check if the mesh is a watertigh triangle mesh
  bool check_mesh(Scene_polyhedron_item* item) {
    double omega_L = ui->omega_L->value();
    double omega_H = ui->omega_H->value();
    double omega_P = ui->omega_P->value();
    double edgelength_TH = ui->edgelength_TH->value();
    double volume_TH = ui->volume_TH->value();
    double area_TH = ui->area_TH->value();
    double diag = scene->len_diagonal();
    bool is_medially_centered = ui->is_medially_centered->isChecked();

    Polyhedron *pMesh = item->polyhedron();

    if (mcs == NULL)
    {
      if (!is_mesh_valid(pMesh))
      {
        return false;
      }

      // save a copy before any operation
      mCopy = new Polyhedron(*pMesh);
      if (is_medially_centered)
      {
        mcs = new Mean_curvature_skeleton(pMesh, Vertex_index_map(), Edge_index_map(),
                                          omega_L, omega_H, omega_P, edgelength_TH, true, volume_TH,
                                          area_TH);
      }
      else
      {
        mcs = new Mean_curvature_skeleton(pMesh, Vertex_index_map(), Edge_index_map(),
                                          omega_L, omega_H, edgelength_TH, volume_TH, area_TH);
      }
      fixedPointsItemIndex = -1;
      nonFixedPointsItemIndex = -1;
      poleLinesItemIndex = -1;

      Scene_polyhedron_item* item_copy = new Scene_polyhedron_item(mCopy);
      copyItemIndex = scene->addItem(item_copy);
      item_copy->setName(QString("original mesh of %1").arg(item->name()));
      item_copy->setVisible(false);
    }
    else
    {
      Polyhedron* mesh = mcs->get_polyhedron();
      if (mesh != pMesh)
      {
        if (!is_mesh_valid(pMesh))
        {
          return false;
        }

        delete mcs;

        // save a copy before any operation
        mCopy = new Polyhedron(*pMesh);
        if (is_medially_centered)
        {
          mcs = new Mean_curvature_skeleton(pMesh, Vertex_index_map(), Edge_index_map(),
                                            omega_L, omega_H, omega_P, edgelength_TH, true, volume_TH,
                                            area_TH);
        }
        else
        {
          mcs = new Mean_curvature_skeleton(pMesh, Vertex_index_map(), Edge_index_map(),
                                            omega_L, omega_H, edgelength_TH, volume_TH, area_TH);
        }
        fixedPointsItemIndex = -1;
        nonFixedPointsItemIndex = -1;
        poleLinesItemIndex = -1;

        Scene_polyhedron_item* item_copy = new Scene_polyhedron_item(mCopy);
        copyItemIndex = scene->addItem(item_copy);
        item_copy->setName(QString("original mesh of %1").arg(item->name()));
        item_copy->setVisible(false);
      }
      else
      {
        mcs->set_omega_L(omega_L);
        mcs->set_omega_H(omega_H);
        mcs->set_omega_P(omega_P);
        mcs->set_edgelength_TH(edgelength_TH);
        mcs->set_volume_TH(volume_TH);
      }
    }
    return true;
  }

public slots:
  void on_actionMCFSkeleton_triggered();
  void on_actionConvert_to_skeleton_triggered();
  void on_actionConvert_to_medial_skeleton_triggered();
  void on_actionContract();
  void on_actionCollapse();
  void on_actionSplit();
  void on_actionDegeneracy();
  void on_actionRun();
  void on_actionSkeletonize();
  void on_actionConverge();
  void on_actionUpdateBBox();
  void on_actionVoronoi();

private:
  Mean_curvature_skeleton* mcs;
  QDockWidget* dockWidget;
  Ui::Mean_curvature_flow_skeleton_plugin* ui;

  int fixedPointsItemIndex;
  int nonFixedPointsItemIndex;
  int poleLinesItemIndex;
  int copyItemIndex;

  Polyhedron *mCopy;
}; // end Polyhedron_demo_mean_curvature_flow_skeleton_plugin

void Polyhedron_demo_mean_curvature_flow_skeleton_plugin::on_actionMCFSkeleton_triggered()
{
  const Scene_interface::Item_id index = scene->mainSelectionIndex();

  Scene_polyhedron_item* item =
    qobject_cast<Scene_polyhedron_item*>(scene->item(index));

  if(item)
  {
    Polyhedron* pMesh = item->polyhedron();

    if(!pMesh) return;

    dockWidget = new QDockWidget(mw);
    ui = new Ui::Mean_curvature_flow_skeleton_plugin();
    ui->setupUi(dockWidget);
    dockWidget->setFeatures(QDockWidget::DockWidgetMovable
                          | QDockWidget::DockWidgetFloatable
                          | QDockWidget::DockWidgetClosable);
    dockWidget->setWindowTitle("Mean Curvature Flow Skeleton");
    mw->addDockWidget(Qt::LeftDockWidgetArea, dockWidget);
    mw->tabifyDockWidget(static_cast<MainWindow*>(mw)->get_ui()->consoleDockWidget, dockWidget);
    dockWidget->show();
    dockWidget->raise();

    connect(ui->pushButton_contract, SIGNAL(clicked()),
            this, SLOT(on_actionContract()));
    connect(ui->pushButton_collapse, SIGNAL(clicked()),
            this, SLOT(on_actionCollapse()));
    connect(ui->pushButton_split, SIGNAL(clicked()),
            this, SLOT(on_actionSplit()));
    connect(ui->pushButton_degeneracy, SIGNAL(clicked()),
            this, SLOT(on_actionDegeneracy()));
    connect(ui->pushButton_run, SIGNAL(clicked()),
            this, SLOT(on_actionRun()));
    connect(ui->pushButton_skeletonize, SIGNAL(clicked()),
            this, SLOT(on_actionSkeletonize()));
    connect(ui->pushButton_converge, SIGNAL(clicked()),
            this, SLOT(on_actionConverge()));
    connect(dynamic_cast<Scene*>(scene), SIGNAL(updated_bbox()),
            this, SLOT(on_actionUpdateBBox()));
    connect(ui->pushButton_voronoi, SIGNAL(clicked()),
            this, SLOT(on_actionVoronoi()));

    double diag = scene->len_diagonal();
    init_ui(diag);

    fixedPointsItemIndex = -1;
    nonFixedPointsItemIndex = -1;
    poleLinesItemIndex = -1;
    copyItemIndex = -1;
  }
}

void Polyhedron_demo_mean_curvature_flow_skeleton_plugin::on_actionUpdateBBox()
{
  double diag = scene->len_diagonal();
  ui->edgelength_TH->setValue(0.002 * diag);
}

void Polyhedron_demo_mean_curvature_flow_skeleton_plugin::on_actionVoronoi()
{
  const Scene_interface::Item_id index = scene->mainSelectionIndex();
  if (!check_item_index(index))
  {
    return;
  }

  Scene_polyhedron_item* item =
    qobject_cast<Scene_polyhedron_item*>(scene->item(index));
  Polyhedron* pMesh = item->polyhedron();

  if (!check_mesh(item))
  {
    return;
  }

  QTime time;
  time.start();
  QApplication::setOverrideCursor(Qt::WaitCursor);

  Scene_polylines_item* poleLinesItem = new Scene_polylines_item();

  std::vector<Point> pole_points;
  mcs->get_poles(pole_points);
  vertex_iterator vb, ve;
  int id = 0;
  for (boost::tie(vb, ve) = boost::vertices(*pMesh); vb != ve; vb++)
  {
    std::vector<Point> line;
    line.clear();

    vertex_descriptor v = *vb;
    Point s = v->point();
    Point t = pole_points[id++];

    line.push_back(s);
    line.push_back(t);
    poleLinesItem->polylines.push_back(line);
  }
  if (poleLinesItemIndex == -1)
  {
    poleLinesItemIndex = scene->addItem(poleLinesItem, false);
  }
  else
  {
    scene->replaceItem(poleLinesItemIndex, poleLinesItem, false);
  }

  scene->itemChanged(index);
  scene->setSelectedItem(index);
  QApplication::restoreOverrideCursor();
}

void Polyhedron_demo_mean_curvature_flow_skeleton_plugin::on_actionConvert_to_skeleton_triggered()
{
  double diag = scene->len_diagonal();
  double omega_L = 1;
  double omega_H = 0.1;
  double edgelength_TH = 0.002 * diag;
  double volume_TH = 1e-04;
  double area_TH = 1e-4;

  const Scene_interface::Item_id index = scene->mainSelectionIndex();

  Scene_polyhedron_item* item =
    qobject_cast<Scene_polyhedron_item*>(scene->item(index));

  if(item)
  {
    Polyhedron* pMesh = item->polyhedron();

    if(!pMesh) return;

    Polyhedron tempMesh = *pMesh;

    Mean_curvature_skeleton* temp_mcs = new Mean_curvature_skeleton(&tempMesh, Vertex_index_map(), Edge_index_map(),
                                      omega_L, omega_H, edgelength_TH, volume_TH, area_TH);

    QTime time;
    time.start();
    QApplication::setOverrideCursor(Qt::WaitCursor);

    temp_mcs->run_to_converge();

    Graph g;
    std::vector<Point> points;

    temp_mcs->convert_to_skeleton();
    temp_mcs->get_skeleton(g, points);

    std::cout << "ok (" << time.elapsed() << " ms, " << ")" << std::endl;

    Scene_polylines_item* skeleton = new Scene_polylines_item();

    boost::graph_traits<Graph>::edge_iterator ei, ei_end;
    for (boost::tie(ei, ei_end) = boost::edges(g); ei != ei_end; ++ei)
    {
      std::vector<Point> line;
      line.clear();
      Point s = points[boost::source(*ei, g)];
      Point t = points[boost::target(*ei, g)];
      line.push_back(s);
      line.push_back(t);
      skeleton->polylines.push_back(line);
    }
    skeleton->setName(QString("skeleton curve of %1").arg(item->name()));
    scene->addItem(skeleton, false);
    item->setGouraudMode();
    item->switch_transparency_on_off();

    QApplication::restoreOverrideCursor();

    delete temp_mcs;
  }
}

void Polyhedron_demo_mean_curvature_flow_skeleton_plugin::on_actionConvert_to_medial_skeleton_triggered()
{
  double diag = scene->len_diagonal();
  double omega_L = 1;
  double omega_H = 0.1;
  double omega_P = 0.2;
  double edgelength_TH = 0.002 * diag;
  double volume_TH = 1e-04;
  double area_TH = 1e-4;

  const Scene_interface::Item_id index = scene->mainSelectionIndex();

  Scene_polyhedron_item* item =
    qobject_cast<Scene_polyhedron_item*>(scene->item(index));

  if(item)
  {
    Polyhedron* pMesh = item->polyhedron();

    if(!pMesh) return;

    Polyhedron tempMesh = *pMesh;

    Mean_curvature_skeleton* temp_mcs = new Mean_curvature_skeleton(&tempMesh, Vertex_index_map(), Edge_index_map(),
                                                                    omega_L, omega_H, omega_P, edgelength_TH, true, volume_TH,
                                                                    area_TH);

    QTime time;
    time.start();
    QApplication::setOverrideCursor(Qt::WaitCursor);

    temp_mcs->run_to_converge();

    Graph g;
    std::vector<Point> points;

    temp_mcs->convert_to_skeleton();
    temp_mcs->get_skeleton(g, points);

    std::cout << "ok (" << time.elapsed() << " ms, " << ")" << std::endl;

    Scene_polylines_item* skeleton = new Scene_polylines_item();

    boost::graph_traits<Graph>::edge_iterator ei, ei_end;
    for (boost::tie(ei, ei_end) = boost::edges(g); ei != ei_end; ++ei)
    {
      std::vector<Point> line;
      line.clear();
      Point s = points[boost::source(*ei, g)];
      Point t = points[boost::target(*ei, g)];
      line.push_back(s);
      line.push_back(t);
      skeleton->polylines.push_back(line);
    }
    skeleton->setName(QString("skeleton curve of %1").arg(item->name()));
    scene->addItem(skeleton, false);
    item->setGouraudMode();
    item->switch_transparency_on_off();

    QApplication::restoreOverrideCursor();

    delete temp_mcs;
  }
}

void Polyhedron_demo_mean_curvature_flow_skeleton_plugin::on_actionContract()
{
  const Scene_interface::Item_id index = scene->mainSelectionIndex();
  if (!check_item_index(index))
  {
    return;
  }

  Scene_polyhedron_item* item =
    qobject_cast<Scene_polyhedron_item*>(scene->item(index));
  Polyhedron* pMesh = item->polyhedron();

  if (!check_mesh(item))
  {
    return;
  }

  QTime time;
  time.start();
  std::cout << "Contract...\n";
  QApplication::setOverrideCursor(Qt::WaitCursor);

  mcs->contract_geometry();

  std::cout << "ok (" << time.elapsed() << " ms, " << ")" << std::endl;

  // update scene
  scene->itemChanged(index);
  QApplication::restoreOverrideCursor();
}

void Polyhedron_demo_mean_curvature_flow_skeleton_plugin::on_actionCollapse()
{
  const Scene_interface::Item_id index = scene->mainSelectionIndex();
  if (!check_item_index(index))
  {
    return;
  }

  Scene_polyhedron_item* item =
    qobject_cast<Scene_polyhedron_item*>(scene->item(index));
  Polyhedron* pMesh = item->polyhedron();

  if (!check_mesh(item))
  {
    return;
  }

  QTime time;
  time.start();
  std::cout << "Collapse...\n";
  QApplication::setOverrideCursor(Qt::WaitCursor);

  std::cout << "before collapse edges\n";
  int num_collapses = mcs->collapse_short_edges();
  std::cout << "collapse " << num_collapses << " edges.\n";

  std::cout << "ok (" << time.elapsed() << " ms, " << ")" << std::endl;

  item->color();
  // update scene
  scene->itemChanged(index);
  QApplication::restoreOverrideCursor();
}

void Polyhedron_demo_mean_curvature_flow_skeleton_plugin::on_actionSplit()
{
  const Scene_interface::Item_id index = scene->mainSelectionIndex();
  if (!check_item_index(index))
  {
    return;
  }

  Scene_polyhedron_item* item =
    qobject_cast<Scene_polyhedron_item*>(scene->item(index));
  Polyhedron* pMesh = item->polyhedron();

  if (!check_mesh(item))
  {
    return;
  }

  QTime time;
  time.start();
  std::cout << "Split...\n";
  QApplication::setOverrideCursor(Qt::WaitCursor);

  std::cout << "before split triangles\n";
  int num_split = mcs->iteratively_split_triangles();
  std::cout << "split " << num_split << " triangles.\n";

  std::cout << "ok (" << time.elapsed() << " ms, " << ")" << std::endl;

  // update scene
  scene->itemChanged(index);
  QApplication::restoreOverrideCursor();
}

void Polyhedron_demo_mean_curvature_flow_skeleton_plugin::on_actionDegeneracy()
{
  const Scene_interface::Item_id index = scene->mainSelectionIndex();
  if (!check_item_index(index))
  {
    return;
  }

  Scene_polyhedron_item* item =
    qobject_cast<Scene_polyhedron_item*>(scene->item(index));
  Polyhedron* pMesh = item->polyhedron();

  if (!check_mesh(item))
  {
    return;
  }

  QTime time;
  time.start();
  std::cout << "Degeneracy\n";
  QApplication::setOverrideCursor(Qt::WaitCursor);

//  mcs->detect_degeneracies_in_disk();
  mcs->detect_degeneracies();

  std::cout << "ok (" << time.elapsed() << " ms, " << ")" << std::endl;

  Scene_points_with_normal_item* fixedPointsItem = new Scene_points_with_normal_item;
  fixedPointsItem->setName(QString("fixed points of %1").arg(item->name()));

  std::vector<Point> fixedPoints;
  mcs->get_fixed_points(fixedPoints);

  Point_set *ps = fixedPointsItem->point_set();
  for (size_t i = 0; i < fixedPoints.size(); i++)
  {
    UI_point_3<Kernel> point(fixedPoints[i].x(), fixedPoints[i].y(), fixedPoints[i].z());
    ps->select(&point);
    ps->push_back(point);
  }

  if (fixedPointsItemIndex == -1)
  {
    fixedPointsItemIndex = scene->addItem(fixedPointsItem, false);
    std::cerr << "add item " << fixedPointsItemIndex << "\n";
  }
  else
  {
    std::cerr << "replace item " << fixedPointsItemIndex << "\n";
    Scene_item* temp = scene->replaceItem(fixedPointsItemIndex, fixedPointsItem, false);
    delete temp;
  }
  // update scene
  scene->itemChanged(index);
  scene->itemChanged(fixedPointsItemIndex);
  scene->setSelectedItem(index);
  QApplication::restoreOverrideCursor();
}

void Polyhedron_demo_mean_curvature_flow_skeleton_plugin::on_actionRun()
{
  const Scene_interface::Item_id index = scene->mainSelectionIndex();
  if (!check_item_index(index))
  {
    return;
  }

  Scene_polyhedron_item* item =
    qobject_cast<Scene_polyhedron_item*>(scene->item(index));
  Polyhedron* pMesh = item->polyhedron();

  if (!check_mesh(item))
  {
    return;
  }

  QTime time;
  time.start();
  QApplication::setOverrideCursor(Qt::WaitCursor);

  mcs->contract();

  std::cout << "ok (" << time.elapsed() << " ms, " << ")" << std::endl;

  // update scene
  Scene_points_with_normal_item* fixedPointsItem = new Scene_points_with_normal_item;
  fixedPointsItem->setName(QString("fixed points of %1").arg(item->name()));

  std::vector<Point> fixedPoints;
  mcs->get_fixed_points(fixedPoints);

  Point_set *ps = fixedPointsItem->point_set();
  for (size_t i = 0; i < fixedPoints.size(); i++)
  {
    UI_point_3<Kernel> point(fixedPoints[i].x(), fixedPoints[i].y(), fixedPoints[i].z());
    ps->select(&point);
    ps->push_back(point);
  }
  if (fixedPointsItemIndex == -1)
  {
    fixedPointsItemIndex = scene->addItem(fixedPointsItem, false);
  }
  else
  {
    Scene_item* temp = scene->replaceItem(fixedPointsItemIndex, fixedPointsItem, false);
    delete temp;
  }

  // draw non-fixed points
//  Scene_points_with_normal_item* nonFixedPointsItem = new Scene_points_with_normal_item;
//  nonFixedPointsItem->setName("non-fixed points");
//  nonFixedPointsItem->setColor(QColor(0, 255, 0));
//  std::vector<Point> nonFixedPoints;
//  mcs->get_non_fixed_points(nonFixedPoints);
//  ps = nonFixedPointsItem->point_set();
//  for (size_t i = 0; i < nonFixedPoints.size(); i++)
//  {
//    UI_point_3<Kernel> point(nonFixedPoints[i].x(), nonFixedPoints[i].y(), nonFixedPoints[i].z());
//    ps->push_back(point);
//  }
//  if (nonFixedPointsItemIndex == -1)
//  {
//    nonFixedPointsItemIndex = scene->addItem(nonFixedPointsItem, false);
//  }
//  else
//  {
//    scene->replaceItem(nonFixedPointsItemIndex, nonFixedPointsItem, false);
//  }

  // draw lines connecting surface points and their correspondent poles
//  Scene_polylines_item* poleLinesItem = new Scene_polylines_item();

//  std::vector<Point> pole_points;
//  mcs->get_poles(pole_points);
//  vertex_iterator vb, ve;
//  int id = 0;
//  for (boost::tie(vb, ve) = boost::vertices(*pMesh); vb != ve; vb++)
//  {
//    std::vector<Point> line;
//    line.clear();

//    vertex_descriptor v = *vb;
//    Point s = v->point();
//    Point t = pole_points[id++];

//    line.push_back(s);
//    line.push_back(t);
//    poleLinesItem->polylines.push_back(line);
//  }
//  if (poleLinesItemIndex == -1)
//  {
//    poleLinesItemIndex = scene->addItem(poleLinesItem, false);
//  }
//  else
//  {
//    scene->replaceItem(poleLinesItemIndex, poleLinesItem, false);
//  }

  scene->itemChanged(index);
  scene->itemChanged(fixedPointsItemIndex);
//  scene->itemChanged(nonFixedPointsItemIndex);
  scene->setSelectedItem(index);
  QApplication::restoreOverrideCursor();
}

void Polyhedron_demo_mean_curvature_flow_skeleton_plugin::on_actionSkeletonize()
{
  const Scene_interface::Item_id index = scene->mainSelectionIndex();
  if (!check_item_index(index))
  {
    return;
  }

  Scene_polyhedron_item* item =
    qobject_cast<Scene_polyhedron_item*>(scene->item(index));
  Polyhedron* pMesh = item->polyhedron();

  if (!check_mesh(item))
  {
    return;
  }

  QTime time;
  time.start();
  QApplication::setOverrideCursor(Qt::WaitCursor);

  Graph g;
  std::vector<Point> points;
  std::vector<std::vector<int> > corr;

  mcs->convert_to_skeleton();
  mcs->get_skeleton(g, points);
  mcs->get_correspondent_vertices(corr);

  std::cout << "ok (" << time.elapsed() << " ms, " << ")" << std::endl;

  Scene_polylines_item* skeleton = new Scene_polylines_item();

  boost::graph_traits<Graph>::edge_iterator ei, ei_end;
  for (boost::tie(ei, ei_end) = boost::edges(g); ei != ei_end; ++ei)
  {
    std::vector<Point> line;
    line.clear();
    Point s = points[boost::source(*ei, g)];
    Point t = points[boost::target(*ei, g)];
    line.push_back(s);
    line.push_back(t);
    skeleton->polylines.push_back(line);
  }
  skeleton->setColor(QColor(255, 0, 0));
  skeleton->setName(QString("skeleton curve of %1").arg(item->name()));
  scene->addItem(skeleton, false);

  vertex_iterator vb, ve;
  std::vector<vertex_descriptor> id_to_vd;
  id_to_vd.clear();
  id_to_vd.resize(boost::num_vertices(*mCopy));
  int id = 0;
  for (boost::tie(vb, ve) = boost::vertices(*mCopy); vb != ve; vb++)
  {
    vertex_descriptor v = *vb;
    id_to_vd[id++] = v;
  }

  Scene_polylines_item* lines = new Scene_polylines_item();

  for (size_t i = 0; i < corr.size(); i++)
  {
    Point s = points[i];
    for (size_t j = 0; j < corr[i].size(); j++)
    {
      std::vector<Point> line;
      line.clear();
      Point t = id_to_vd[corr[i][j]]->point();
      line.push_back(s);
      line.push_back(t);
      lines->polylines.push_back(line);
    }
  }
  lines->setName(QString("correspondent vertices of %1").arg(item->name()));
  lines->setVisible(false);
  scene->addItem(lines, false);

  // set the fixed points and contracted mesh as invisible
  if (fixedPointsItemIndex >= 0)
  {
    scene->item(fixedPointsItemIndex)->setVisible(false);
  }
  // display the original mesh in transparent mode
  item->setVisible(false);
  if (copyItemIndex >= 0)
  {
    scene->item(copyItemIndex)->setVisible(true);
    dynamic_cast<Scene_polyhedron_item*>(scene->item(copyItemIndex))->switch_transparency_on_off();
    scene->item(copyItemIndex)->setGouraudMode();
  }

  // display the end points and junction points
  Scene_points_with_normal_item* endPointsItem = new Scene_points_with_normal_item;
  endPointsItem->setName(QString("end points of %1").arg(item->name()));
  Scene_points_with_normal_item* junctionPointsItem = new Scene_points_with_normal_item;
  junctionPointsItem->setName(QString("junction points of %1").arg(item->name()));

  Point_set *end_ps = endPointsItem->point_set();
  end_ps->set_selected_color(QColor(0, 0, 255));
  end_ps->set_selected_diameter(6.0);
  Point_set *junction_ps = junctionPointsItem->point_set();
  junction_ps->set_selected_color(QColor(0, 255, 0));
  junction_ps->set_selected_diameter(6.0);

  boost::graph_traits<Graph>::vertex_iterator vi;
  for (vi = vertices(g).first; vi != vertices(g).second; ++vi)
  {
    int deg = boost::out_degree(*vi, g);
    if (deg == 1)
    {
      UI_point_3<Kernel> point(points[*vi].x(), points[*vi].y(), points[*vi].z());
      end_ps->select(&point);
      end_ps->push_back(point);
    }
    else if (deg > 2)
    {
      UI_point_3<Kernel> point(points[*vi].x(), points[*vi].y(), points[*vi].z());
      junction_ps->select(&point);
      junction_ps->push_back(point);
    }
  }

  scene->addItem(endPointsItem, false);
  scene->addItem(junctionPointsItem, false);

  // add segmentation
  std::vector<int> skeleton_segment;
  std::vector<bool> deleted;
  skeleton_segment.resize(boost::num_vertices(g));
  deleted.resize(boost::num_vertices(g), false);
  for (vi = vertices(g).first; vi != vertices(g).second; ++vi)
  {
    int deg = boost::out_degree(*vi, g);
    if (deg > 2)
    {
      // for branching point, cut some incident edges to make it an end point
      out_edge_iter e, e_end;
      deleted[*vi] = true;
      bool move_corr = false;
      for (boost::tie(e, e_end) = boost::out_edges(*vi, g); e != e_end; e++)
      {
        edge_desc ed = *e;
        int target = boost::target(ed, g);
        // delete the branching point and move correspondent vertices to another vertex
        if (!move_corr && !deleted[target])
        {
          corr[target].insert(corr[target].end(),
                              corr[*vi].begin(),
                              corr[*vi].end());
          move_corr = true;
          break;
        }
      }
    }
  }

  int num_segment = 0;
  std::vector<bool> visited;
  visited.resize(boost::num_vertices(g), false);
  for (vi = vertices(g).first; vi != vertices(g).second; ++vi)
  {
    int vid = *vi;
    std::queue<int> qu;
    while(!qu.empty())
    {
      qu.pop();
    }

    // branching points have been deleted
    if (deleted[vid])
    {
      continue;
    }

    if (!visited[vid])
    {
      qu.push(vid);
      while (!qu.empty())
      {
        int cur = qu.front();
        qu.pop();
        visited[cur] = true;
        skeleton_segment[cur] = num_segment;
        out_edge_iter e, e_end;
        for (boost::tie(e, e_end) = boost::out_edges(cur, g); e != e_end; e++)
        {
          edge_desc ed = *e;
          int target = boost::target(ed, g);
          if (!visited[target] && !deleted[target])
          {
            qu.push(target);
          }
        }
      }
      num_segment++;
    }
  }

  std::cout << "num segment " << num_segment << "\n";

  std::vector<int> segment_id;
  segment_id.resize(boost::num_vertices(*mCopy), -1);
  for (vi = vertices(g).first; vi != vertices(g).second; ++vi)
  {
    int vid = *vi;
    int seg = skeleton_segment[vid];
    for (size_t i = 0; i < corr[vid].size(); i++)
    {
      segment_id[corr[vid][i]] = seg;
    }
  }

  Polyhedron *segment_mesh = new Polyhedron(*mCopy);
  int vertex_id_count = 0;
  for (boost::tie(vb, ve) = boost::vertices(*segment_mesh); vb != ve; ++vb)
  {
    vb->id() = vertex_id_count++;
  }

  for (Facet_iterator f = segment_mesh->facets_begin(); f != segment_mesh->facets_end(); ++f)
  {
    Polyhedron::Halfedge_const_handle he = f->facet_begin();
    int vid1 = he->vertex()->id();
    int vid2 = he->next()->vertex()->id();
    int vid3 = he->next()->next()->vertex()->id();
    int sid1 = segment_id[vid1];
    int sid2 = segment_id[vid2];
    int sid3 = segment_id[vid3];

    int id;
    if (sid1 == sid2 && sid2 == sid3)
    {
      id = sid1;
    }
    else if (sid1 == sid2)
    {
      id = sid1;
    }
    else if (sid1 == sid3)
    {
      id = sid1;
    }
    else if (sid2 == sid3)
    {
      id = sid3;
    }
    else
    {
      id = sid1;
    }

    f->set_patch_id(id);
  }

  Scene_polyhedron_item* item_segmentation = new Scene_polyhedron_item(segment_mesh);
  scene->addItem(item_segmentation);
  item_segmentation->setName(QString("segmentation of %1").arg(item->name()));
  item_segmentation->setVisible(false);

  // update scene
  QApplication::restoreOverrideCursor();
}

void Polyhedron_demo_mean_curvature_flow_skeleton_plugin::on_actionConverge()
{
  const Scene_interface::Item_id index = scene->mainSelectionIndex();
  if (!check_item_index(index))
  {
    return;
  }

  Scene_polyhedron_item* item =
    qobject_cast<Scene_polyhedron_item*>(scene->item(index));
  Polyhedron* pMesh = item->polyhedron();

  if (!check_mesh(item))
  {
    return;
  }

  QTime time;
  time.start();
  QApplication::setOverrideCursor(Qt::WaitCursor);

  mcs->run_to_converge();

  std::cout << "ok (" << time.elapsed() << " ms, " << ")" << std::endl;

  // update scene
  Scene_points_with_normal_item* fixedPointsItem = new Scene_points_with_normal_item;
  fixedPointsItem->setName(QString("fixed points of %1").arg(item->name()));

  std::vector<Point> fixedPoints;
  mcs->get_fixed_points(fixedPoints);

  Point_set *ps = fixedPointsItem->point_set();
  for (size_t i = 0; i < fixedPoints.size(); i++)
  {
    UI_point_3<Kernel> point(fixedPoints[i].x(), fixedPoints[i].y(), fixedPoints[i].z());
    ps->select(&point);
    ps->push_back(point);
  }
  if (fixedPointsItemIndex == -1)
  {
    fixedPointsItemIndex = scene->addItem(fixedPointsItem, false);
  }
  else
  {
    Scene_item* temp = scene->replaceItem(fixedPointsItemIndex, fixedPointsItem, false);
    delete temp;
  }

//  Scene_points_with_normal_item* nonFixedPointsItem = new Scene_points_with_normal_item;
//  nonFixedPointsItem->setName("non-fixed points");
//  nonFixedPointsItem->setColor(QColor(0, 255, 0));
//  std::vector<Point> nonFixedPoints;
//  mcs->get_non_fixed_points(nonFixedPoints);
//  ps = nonFixedPointsItem->point_set();
//  for (size_t i = 0; i < nonFixedPoints.size(); i++)
//  {
//    UI_point_3<Kernel> point(nonFixedPoints[i].x(), nonFixedPoints[i].y(), nonFixedPoints[i].z());
//    ps->push_back(point);
//  }
//  if (nonFixedPointsItemIndex == -1)
//  {
//    nonFixedPointsItemIndex = scene->addItem(nonFixedPointsItem, false);
//  }
//  else
//  {
//    scene->replaceItem(nonFixedPointsItemIndex, nonFixedPointsItem, false);
//  }

  scene->itemChanged(index);
  scene->itemChanged(fixedPointsItemIndex);
//  scene->itemChanged(nonFixedPointsItemIndex);
  scene->setSelectedItem(index);
  QApplication::restoreOverrideCursor();
}

Q_EXPORT_PLUGIN2(Polyhedron_demo_mean_curvature_flow_skeleton_plugin, Polyhedron_demo_mean_curvature_flow_skeleton_plugin)

#include "Polyhedron_demo_mean_curvature_flow_skeleton_plugin.moc"
