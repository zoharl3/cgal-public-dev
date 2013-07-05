#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>

#include <cassert>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

#define SSTR( x ) dynamic_cast< std::ostringstream & >( \
        ( std::ostringstream() << std::dec << x ) ).str()



typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

typedef CGAL::Triangulation_vertex_base_2<K>                     Vb;
typedef CGAL::Constrained_triangulation_face_base_2<K>           Fb;
typedef CGAL::Triangulation_data_structure_2<Vb,Fb>              TDS;
typedef CGAL::Exact_predicates_tag                               Itag;
typedef CGAL::Constrained_Delaunay_triangulation_2<K, TDS, Itag> CDT;
typedef CDT::Point          Point;

typedef CDT::Finite_faces_iterator    Finite_faces_iterator;
typedef CDT::Polygon_2                    Polygon_2;

void Triangulation_to_vtk(CDT cdt,std::string name)
{
    Finite_faces_iterator ffi2=cdt.finite_faces_begin();
    std::ofstream out(name.c_str());
    out<<"# vtk DataFile Version 2.0 con "<<cdt.number_of_faces()<<" triangulos y "<<cdt.number_of_vertices()<<" vertices"<<std::endl;
    out<<"grid description"<<std::endl;
    out<<"ASCII"<<std::endl;
    out<<"DATASET UNSTRUCTURED_GRID"<<std::endl<<std::endl;
    out<<"POINTS "<<3*cdt.number_of_faces()<<" double"<<std::endl;

    for(; ffi2!=cdt.finite_faces_end(); ffi2++)
    {
        out<<ffi2->vertex(0)->point().x()<<".0 "<<ffi2->vertex(0)->point().y()<<".0 0.0"<<std::endl;
        out<<ffi2->vertex(1)->point().x()<<".0 "<<ffi2->vertex(1)->point().y()<<".0 0.0"<<std::endl;
        out<<ffi2->vertex(2)->point().x()<<".0 "<<ffi2->vertex(2)->point().y()<<".0 0.0"<<std::endl;
    }
    out<<std::endl<<"CELLS "<<cdt.number_of_faces()<<" "<<cdt.number_of_faces()*4<<std::endl;
    for(int z=0; z<cdt.number_of_faces();z++)
        out<<3<<" "<<z*3<<" "<<z*3+1<<" "<<z*3+2<<std::endl;
    out<<std::endl<<"CELL_TYPES "<<cdt.number_of_faces()<<std::endl;
    for(int g=0;g<cdt.number_of_faces() - 1;g++)
        out<<5<<" ";
    out<<5<<std::endl<<std::endl;
    out<<"CELL_DATA "<<cdt.number_of_faces()<<std::endl;
    out<<"SCALARS scalars float 1"<<std::endl;
    out<<"LOOKUP_TABLE lut"<<std::endl;
    ffi2=cdt.finite_faces_begin();
    
    for(; ffi2!=cdt.finite_faces_end(); ffi2++)
    {
            out<<"1"<<std::endl;
    }
    
    out<<std::endl<<"LOOKUP_TABLE lut 101"<<std::endl;
    
    for(float i=0.00;i<=1.00;i=i+0.01)
    {
        out<<i<<" "<<i<<" "<<i<<" 1.0"<<std::endl;
    }
}


void voronoi_cells_to_vtk(Polygon_2 poly, std::string name){
  if (poly.size() == 0){
    return;
  }
  std::ofstream out(name.c_str());
  out<<"# vtk DataFile Version 1.0"<<std::endl;
  out<<"Line representation of vtk"<<std::endl;
  out<<"ASCII"<<std::endl<<std::endl;
  out<<"DATASET POLYDATA"<<std::endl;
  out<<"POINTS "<<poly.size()<<" float"<<std::endl;
  for(int i=0;i<poly.size();i++){
    out<<poly[i]<<" 0.0"<<std::endl;
  }
  out<<std::endl;
  out<<"LINES "<<"1 "<<poly.size()+2<<std::endl<<poly.size()+1;
  for(int i=1;i<=poly.size();i++){
    out<<" "<<i;
  }
  out<<" 1\n";
}

int main(int argc, char **argv)
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

  std::cout << "The number of vertices is  ";
  std::cout << cdt.number_of_vertices() << std::endl;

  std::cout << "Build the constrained Voronoi diagram" << std::endl;
  std::cout << "Output CVD to an image" << std::endl;
  Triangulation_to_vtk(cdt,"./vtk_files/cdt.vtk");

  int i=1;
  for(CDT::Finite_vertices_iterator vit = cdt.finite_vertices_begin();
      vit != cdt.finite_vertices_end();
      ++vit)
  {
    if (!cdt.cell_is_infinite(vit)){
      Polygon_2 poly = cdt.dual(vit);
      std::string name = "./vtk_files/voronoi_cell";
      name.append(SSTR(i));
      name.append(".vtk");
      voronoi_cells_to_vtk(poly,name);
      i++;
    }
    
  }

  return 0;
}
