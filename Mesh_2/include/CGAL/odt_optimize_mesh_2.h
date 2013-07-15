#ifndef CGAL_ODT_OPTIMIZE_MESH_2_H
#define CGAL_ODT_OPTIMIZE_MESH_2_H


#include <CGAL/Mesh_2/Mesh_global_optimizer.h>
#include <CGAL/Mesh_2/Odt_move.h>

namespace CGAL
{
  template <typename CDT> 
  void
  odt_optimize_mesh_2(CDT& cdt,
                      const std::size_t& max_iteration_number)
{
  typedef typename Mesh_2::Odt_move<CDT>    Move;
  typedef typename Mesh_2::Mesh_global_optimizer<CDT,Move> Odt_optimizer;
  
  // Create optimizer
  Odt_optimizer opt(cdt);
   
  // 1000 iteration max to avoid infinite loops
  if ( 0 == max_iteration_number )
    max_iteration_number = 1000;
  
  // Launch optimization
  return opt(static_cast<int>(max_iteration_number));
}
}//namespace CGAL

#endif