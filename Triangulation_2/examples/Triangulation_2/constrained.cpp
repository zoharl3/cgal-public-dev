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

int
main( )
{
  CDT cdt;
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

/*  CDT cdt;
  cdt.insert(Point(0,0));
  cdt.insert(Point(0,5));
  cdt.insert(Point(5,0));
  cdt.insert(Point(5,5));

  for(CDT::Finite_faces_iterator ffi = cdt.finite_faces_begin();
        ffi != cdt.finite_faces_end(); ffi++){
    Edge edge0 = cdt.mirror_edge(Edge(ffi,0));
    //std::cout << "is infinite " << cdt.is_infinite((edge0.first)->vertex(edge0.second)) << std::endl;
    std::cout << "is constrained " << cdt.is_constrained(Edge(ffi,0)) << std::endl;
    Edge edge1 = cdt.mirror_edge(Edge(ffi,1));
    //std::cout << "is infinite " << cdt.is_infinite((edge1.first)->vertex(edge1.second)) << std::endl;
    std::cout << "is constrained " << cdt.is_constrained(Edge(ffi,1)) << std::endl;
    Edge edge2 = cdt.mirror_edge(Edge(ffi,2));
    //std::cout << "is infinite " << cdt.is_infinite((edge2.first)->vertex(edge2.second)) << std::endl;
    std::cout << "is constrained " << cdt.is_constrained(Edge(ffi,2)) << std::endl;
  }*/
  return 0;
}
