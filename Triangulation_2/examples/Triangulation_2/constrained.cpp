#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>

#include <cassert>
#include <iostream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

typedef CGAL::Triangulation_vertex_base_2<K>                     Vb;
typedef CGAL::Constrained_triangulation_face_base_2<K>           Fb;
typedef CGAL::Triangulation_data_structure_2<Vb,Fb>              TDS;
typedef CGAL::Exact_predicates_tag                               Itag;
typedef CGAL::Constrained_Delaunay_triangulation_2<K, TDS, Itag> CDT;
typedef CDT::Point          Point;
typedef CDT::Edge           Edge;
typedef CDT::Vertex_handle  Vertex_handle;
typedef CDT::Finite_vertices_iterator Finite_vertices_iterator;
typedef CDT::Polygon_2  Polygon;

int
main( )
{
/*  CDT cdt;
  std::cout << "Inserting a grid of 5x5 constraints " << std::endl;
  for (int i = 1; i < 6; ++i)
    cdt.insert_constraint( Point(0,i), Point(6,i));
  for (int j = 1; j < 6; ++j)
    cdt.insert_constraint( Point(j,0), Point(j,6));

  assert(cdt.is_valid());
  int count = 0;
  for (CDT::Finite_edges_iterator eit = cdt.finite_edges_begin();
       eit != cdt.finite_edges_end();
       ++eit)
    if (cdt.is_constrained(*eit)) ++count;
  std::cout << "The number of resulting constrained edges is  ";
  std::cout <<  count << std::endl;
*/
  CDT cdt;
  //cdt.insert(Point(0,0));
  cdt.insert(Point(0,5));
  cdt.insert(Point(5,0));
  cdt.insert(Point(-5,0));
  cdt.insert(Point(0,-5));
  Vertex_handle vh = cdt.insert(Point(0,0));
  Vertex_handle vh2 = cdt.insert(Point(5,5));
  //cdt.remove(vh);
  for(Finite_vertices_iterator vit = cdt.finite_vertices_begin();
    vit != cdt.finite_vertices_end();
    vit++){
    Polygon poly = cdt.dual(vit);
    std::cout<<poly<<std::endl;  
  }
  //std::cout<<"0 = "<<cdt.cell_is_infinite(vh)<<std::endl;
  std::cout<<"0 = "<<cdt.cell_is_infinite(vh2)<<std::endl;
  return 0;
}
