#include <CGAL/Hole_filling.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <vector>
#include <boost/tuple/tuple.hpp>
#include <CGAL/utility.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3 Point_3;

struct My_triangle {
  int v0, v1, v2;
  My_triangle(int v0, int v1, int v2) : v0(v0), v1(v1), v2(v2) { }
};

int main() {
  std::vector<Point_3> polyline;
  polyline.push_back(Point_3( 1.,0.,0.));
  polyline.push_back(Point_3( 0.,1.,0.));
  polyline.push_back(Point_3(-1.,0.,0.));
  polyline.push_back(Point_3( 1.,1.,0.));
  // repeating first point (i.e. polyline.push_back(Point_3(1.,0.,0.)) ) is optional

  // Types having (int, int, int) constructor can be used to hold output triangles
  std::vector<boost::tuple<int, int, int> > patch_1;
  std::vector<CGAL::Triple<int, int, int> > patch_2;
  std::vector<My_triangle>                  patch_3;

  patch_1.reserve(polyline.size() -2); // there will be exactly n-2 triangles in the patch
  CGAL::triangulate_hole_polyline(polyline.begin(), polyline.end(), back_inserter(patch_1));
  CGAL::triangulate_hole_polyline(polyline.begin(), polyline.end(), back_inserter(patch_2));
  CGAL::triangulate_hole_polyline(polyline.begin(), polyline.end(), back_inserter(patch_3));

  for(std::size_t i = 0; i < patch_1.size(); ++i) {
    std::cout << "Triangle " << i << ": " << patch_1[i].get<0>() << " " 
              << patch_1[i].get<1>() << " " << patch_1[i].get<2>() << std::endl;

    assert(patch_1[i].get<0>() == patch_2[i].first  && patch_2[i].first  == patch_3[i].v0);
    assert(patch_1[i].get<1>() == patch_2[i].second && patch_2[i].second == patch_3[i].v1);
    assert(patch_1[i].get<2>() == patch_2[i].third  && patch_2[i].third  == patch_3[i].v2);
  }
}
