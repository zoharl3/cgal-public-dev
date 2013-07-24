#define CGAL_MESH_2_OPTIMIZER_VERBOSE

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Delaunay_mesher_2.h>
#include <CGAL/Delaunay_mesh_face_base_2.h>
#include <CGAL/Delaunay_mesh_size_criteria_2.h>

#include <iostream>
#include <fstream>

#define SSTR( x ) dynamic_cast< std::ostringstream & >( \
        ( std::ostringstream() << std::dec << x ) ).str()


typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Triangulation_vertex_base_2<K> Vb;
typedef CGAL::Delaunay_mesh_face_base_2<K> Fb;
typedef CGAL::Triangulation_data_structure_2<Vb, Fb> Tds;
typedef CGAL::Constrained_Delaunay_triangulation_2<K, Tds> CDT;
typedef CGAL::Delaunay_mesh_size_criteria_2<CDT> Criteria;
typedef CGAL::Delaunay_mesher_2<CDT, Criteria> Mesher;
typedef CDT::Finite_faces_iterator                              Finite_faces_iterator;

typedef CDT::Vertex_handle Vertex_handle;
typedef CDT::Point Point;


void Triangulation_to_vtk(const CDT& cdt,const std::string& name)
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
        out<<ffi2->vertex(0)->point().x()<<" "<<ffi2->vertex(0)->point().y()<<" 0.0"<<std::endl;
        out<<ffi2->vertex(1)->point().x()<<" "<<ffi2->vertex(1)->point().y()<<" 0.0"<<std::endl;
        out<<ffi2->vertex(2)->point().x()<<" "<<ffi2->vertex(2)->point().y()<<" 0.0"<<std::endl;
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


int main()
{
  CDT cdt;

  Vertex_handle va = cdt.insert(Point(-4,0));
  Vertex_handle vb = cdt.insert(Point(0,-1));
  Vertex_handle vc = cdt.insert(Point(4,0));
  Vertex_handle vd = cdt.insert(Point(0,1));
  cdt.insert(Point(2, 0.6));

  cdt.insert_constraint(va, vb);
  cdt.insert_constraint(vb, vc);
  cdt.insert_constraint(vc, vd);
  cdt.insert_constraint(vd, va);

  std::cout << "Number of vertices: " << cdt.number_of_vertices() << std::endl;
  Triangulation_to_vtk(cdt,"odt_before_mesher.vtk");

  Mesher mesher(cdt);
  std::cout << "Meshing..." << std::endl;
  // 0.125 is the default shape bound. It corresponds to abound 20.6 degree.
  // 0.5 is the upper bound on the length of the longuest edge.
  // See reference manual for Delaunay_mesh_size_traits_2<K>.
  mesher.set_criteria(Criteria(0.125, 0.5));
  mesher.refine_mesh();

  std::cout << "Number of vertices: " << cdt.number_of_vertices() << std::endl;
  Triangulation_to_vtk(cdt,"odt_mid_mesher.vtk");

  int nb_iterations = 100;
  mesher.odt(nb_iterations);
  std::cout << "Number of vertices: " << cdt.number_of_vertices() << std::endl;
  Triangulation_to_vtk(cdt,"odt_after_mesher.vtk");
}
