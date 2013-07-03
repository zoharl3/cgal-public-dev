#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Projection_traits_xy_3.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Polygon_2.h>

#include <fstream>
#include <vector>


typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Delaunay_triangulation_2<K> Delaunay;

typedef Delaunay::Point                                             Point;
typedef CGAL::Polygon_2<K,std::vector<Point> > Polygon;

int main()
{
  std::ifstream in("data/terrain.cin");
  std::istream_iterator<Point> begin(in);
  std::istream_iterator<Point> end;


  Delaunay dt(begin, end);
  std::cout << dt.number_of_vertices() << std::endl;
  

  Delaunay::Finite_faces_iterator fit= dt.finite_faces_begin();
  std::cout << dt.dual(fit) << std::endl;
  Delaunay::Finite_vertices_iterator vit=dt.finite_vertices_begin();
  vit++;
  Delaunay::Face_circulator fc = dt.incident_faces(vit),done_f(fc);
  do
  {

    if(!dt.is_infinite(fc)){
      std::cout<<"circumcenter: "<<dt.dual(fc)<<std::endl;
    }else{
      std::cout<<"sorpresa"<<std::endl;
    }
  } while (++fc != done_f);

  Polygon poly = dt.dual(vit);

  std::cout << "///////////// Polygon ////////////"<< std::endl;
  std::cout << poly << std::endl;


  return 0;
}
