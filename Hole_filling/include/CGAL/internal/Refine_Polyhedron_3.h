#ifndef CGAL_HOLE_FILLING_REFINE_POLYHEDRON_3_H
#define CGAL_HOLE_FILLING_REFINE_POLYHEDRON_3_H
#include <cmath>
#include <map>
#include <set>
#include <CGAL/assertions.h>
#include <CGAL/trace.h>
#include <CGAL/squared_distance_3.h>

namespace CGAL {
namespace internal {

template<class Polyhedron>
class Refine_Polyhedron_3 {
// typedefs
  typedef typename Polyhedron::Traits::Point_3 Point_3;
  typedef typename Polyhedron::Vertex_handle Vertex_handle;
  typedef typename Polyhedron::Halfedge_handle Halfedge_handle;
  typedef typename Polyhedron::Facet_handle Facet_handle;
  typedef typename Polyhedron::Vertex_iterator Vertex_iterator;
  typedef typename Polyhedron::Halfedge_iterator Halfedge_iterator;
  typedef typename Polyhedron::Halfedge_around_facet_circulator  Halfedge_around_facet_circulator;
  typedef typename Polyhedron::Halfedge_around_vertex_circulator  Halfedge_around_vertex_circulator;

private:
  bool relax(Polyhedron& poly, Halfedge_handle h)
  {
    const Point_3& p = h->vertex()->point();
    const Point_3& q = h->opposite()->vertex()->point();
    const Point_3& r = h->next()->vertex()->point();
    const Point_3& s = h->opposite()->next()->vertex()->point();
    if( (CGAL::ON_UNBOUNDED_SIDE  != CGAL::side_of_bounded_sphere(p,q,r,s)) ||
      (CGAL::ON_UNBOUNDED_SIDE  != CGAL::side_of_bounded_sphere(p,q,s,r)) ){

      //if(!poly.is_valid()) { 
      //  std::cout << "before flip not valid" << std::endl; 
      //}

      // this is the part which makes is_valid not valid
      for(Halfedge_iterator hb = poly.halfedges_begin(); hb != poly.halfedges_end(); ++hb) {
        if(hb->vertex() == hb->next()->vertex())
        {
          std::cout << "before flip not valid" << std::endl;
        }
      }

      poly.flip_edge(h);

      //if(!poly.is_valid()) {
      //  std::cout << "after flip not valid" << std::endl; 
      //}

      // this is the part which makes is_valid not valid
      for(Halfedge_iterator hb = poly.halfedges_begin(); hb != poly.halfedges_end(); ++hb) {
        if(hb->vertex() == hb->next()->vertex())
        {
          std::cout << "after flip not valid" << std::endl;
        }
      }


      return true;
    }
    return false;
  }

  template<class VertexOutputIterator>
  bool subdivide(Polyhedron& poly, 
    std::vector<Facet_handle>& facets, 
    std::set<Facet_handle>& interior_map,
    std::map<Vertex_handle, double>& scale_attribute, 
    VertexOutputIterator vertex_out,
    double alpha)
  {
    std::list<Facet_handle> new_facets;
    for(typename std::vector<Facet_handle>::iterator it = facets.begin(); it!= facets.end(); ++it){
      CGAL_assertion(*it  != Facet_handle());

      Halfedge_handle hh =  (*it)->halfedge();
      Vertex_handle vi = (*it)->halfedge()->vertex();
      Vertex_handle vj = (*it)->halfedge()->next()->vertex();
      Vertex_handle vk = (*it)->halfedge()->prev()->vertex();
      Point_3 c = CGAL::centroid(vi->point(), vj->point(), vk->point());
      double sac  = (scale_attribute[vi] + scale_attribute[vj] + scale_attribute[vk])/3.0;
      double dist_c_vi = std::sqrt(CGAL::squared_distance(c,vi->point()));
      double dist_c_vj = std::sqrt(CGAL::squared_distance(c,vj->point()));
      double dist_c_vk = std::sqrt(CGAL::squared_distance(c,vk->point()));
      if((alpha * dist_c_vi > sac) &&
        (alpha * dist_c_vj > sac) &&
        (alpha * dist_c_vk > sac) &&
        (alpha * dist_c_vi > scale_attribute[vi]) &&
        (alpha * dist_c_vj > scale_attribute[vj]) &&
        (alpha * dist_c_vk > scale_attribute[vk])){
          Halfedge_handle h = poly.create_center_vertex((*it)->halfedge());
          h->vertex()->point() = c;
          scale_attribute[h->vertex()] = sac;
          *vertex_out++ = h->vertex();

          // collect 2 new facets for next round 
          Facet_handle h1 = h->next()->opposite()->face();
          Facet_handle h2 = h->opposite()->face();
          new_facets.push_back(h1); interior_map.insert(h1);
          new_facets.push_back(h2); interior_map.insert(h2);
          // relax edges of the  patching mesh 
          Halfedge_handle e_ij = h->prev();
          Halfedge_handle e_ik = h->opposite()->next();
          Halfedge_handle e_jk = h->next()->opposite()->prev();

          if(interior_map.find(e_ij->opposite()->face()) != interior_map.end()){
            relax(poly, e_ij);
          }
          if(interior_map.find(e_ik->opposite()->face()) != interior_map.end()){
            relax(poly, e_ik);
          }
          if(interior_map.find(e_jk->opposite()->face()) != interior_map.end()){
            relax(poly, e_jk);
          }
      }
    }
    facets.insert(facets.end(), new_facets.begin(), new_facets.end());
    return ! new_facets.empty();
  }

  bool relax(Polyhedron& poly, const std::vector<Facet_handle>& facets, const std::set<Facet_handle>& interior_map)
  {
    int flips = 0;
    std::list<Halfedge_handle> interior_edges;
    std::set<Halfedge_handle> included_map; // do not use just std::set, the order effects the output (for the same input we want to get same output)

    for(typename std::vector<Facet_handle>::const_iterator it = facets.begin(); it!= facets.end(); ++it){
      Halfedge_around_facet_circulator  circ = (*it)->facet_begin(), done(circ);
      do {
        Halfedge_handle h = circ;
        Halfedge_handle oh = h->opposite();
        if(interior_map.find(oh->face()) != interior_map.end()){
          // it's an interior edge
          Halfedge_handle h_rep = (h < oh) ? h : oh;
          if(included_map.insert(h_rep).second) {
            interior_edges.push_back(h_rep);
          }
        }
        ++circ;
      } while(circ != done);
    }

    CGAL_TRACE_STREAM << "Test " << interior_edges.size() << " edges " << std::endl;
    for(typename std::list<Halfedge_handle>::iterator it = interior_edges.begin();
      it != interior_edges.end();
      ++it){
        if(relax(poly,*it)){
          ++flips;
        }
    }

    CGAL_TRACE_STREAM << "|flips| = " << flips << std::endl;
    return flips > 0;
  }

  double average_length(Vertex_handle vh, const std::set<Facet_handle>& interior_map)
  {
    const Point_3& vp = vh->point(); 
    Halfedge_around_vertex_circulator circ(vh->vertex_begin()), done(circ);
    int deg = 0;
    double sum = 0;
    do {
      Facet_handle f(circ->facet()), f_op(circ->opposite()->facet());
      if(interior_map.find(f) != interior_map.end() && interior_map.find(f_op) != interior_map.end())
      { continue; } // which means current edge is an interior edge and should not be included in scale attribute calculation

      const Point_3& vq = circ->opposite()->vertex()->point();
      sum += std::sqrt(CGAL::squared_distance(vp, vq));
      ++deg;
    } while(++circ != done);

    CGAL_assertion(deg != 0); // interior vertices are not accepted
    return sum/deg;
  }

  void calculate_scale_attribute(
    const std::vector<Facet_handle>& facets, 
    const std::set<Facet_handle>& interior_map,
    std::map<Vertex_handle, double>& scale_attribute) 
  {
    for(std::vector<Facet_handle>::const_iterator f_it = facets.begin(); f_it != facets.end(); ++f_it) {
      Halfedge_around_facet_circulator circ((*f_it)->facet_begin()), done(circ);
      do {
        Vertex_handle v = circ->vertex();
        std::pair<std::map<Vertex_handle, double>::iterator, bool> v_insert 
          = scale_attribute.insert(std::make_pair(v, 0));
        if(!v_insert.second) { continue; } // already calculated
        v_insert.first->second = average_length(v, interior_map);
      } while(++circ != done);
    }
  }
public:
  template<class InputIterator, class FacetOutputIterator, class VertexOutputIterator>
  void refine(Polyhedron& poly, 
    InputIterator facet_begin, 
    InputIterator facet_end, 
    FacetOutputIterator facet_out,
    VertexOutputIterator vertex_out,
    double alpha)
  {
    
    std::vector<Facet_handle> facets(facet_begin, facet_end); // do not use just std::set, the order effects the output (for the same input we want to get same output)
    std::set<Facet_handle> interior_map(facet_begin, facet_end);

    std::map<Vertex_handle, double> scale_attribute;
    calculate_scale_attribute(facets, interior_map, scale_attribute);

    int i = 0;
    do {
      if(i == 10){
        break;
      }
      i++;
    } while( subdivide(poly, facets, interior_map, scale_attribute, vertex_out, alpha) && relax(poly,facets, interior_map) );

    std::copy(facets.begin(), facets.end(), facet_out);
    // according to paper it should be like below (?) IOY
    //while(true) {
    //  bool subdiv = subdivide(poly, facets);
    //  if(!subdiv) { break; }
    //  while(relax(poly,facets)) {}
    //}
  }
};

}//namespace internal

/** 
 * @brief Function refining a region on surface mesh.
 *
 * @tparam Polyhedron a %CGAL polyhedron
 * @tparam InputIterator iterator over input facets
 * @tparam FacetOutputIterator iterator holding 'Polyhedron::Facet_handle' for patch facets.
 * @tparam VertexOutputIterator iterator holding 'Polyhedron::Vertex_handle' for patch vertices.
 *
 * @param polyhedron surface mesh to be refined
 * @param facet_begin first iterator of the range of facets
 * @param facet_end past-the-end iterator of the range of facets
 * @param[out] facet_out iterator over patch facets
 * @param[out] vertex_out iterator over patch vertices without including boundary
 * @param density_control_factor factor for density where larger values cause denser refinements
 */
template<class Polyhedron, class InputIterator, class FacetOutputIterator, class VertexOutputIterator>
void refine(Polyhedron& poly, 
  InputIterator facet_begin, 
  InputIterator facet_end,
  FacetOutputIterator facet_out,
  VertexOutputIterator vertex_out,
  double density_control_factor = std::sqrt(2.0)
  )
{
  internal::Refine_Polyhedron_3<Polyhedron> refine_functor;
  refine_functor.refine(poly, facet_begin, facet_end, facet_out, vertex_out, density_control_factor);
}



}//namespace CGAL
#endif //CGAL_HOLE_FILLING_REFINE_POLYHEDRON_3_H