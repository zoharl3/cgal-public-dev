namespace CGAL {

/*!
\ingroup PkgMesh2Functions

The function `lloyd_optimize_mesh_2()` is a mesh optimization process 
based on the minimization of a global energy function. 

In `lloyd_optimize_mesh_2()`, the minimized global energy may be interpreted 
as the \f$ L^1\f$-norm of the error achieved 
when the function \f$ x^2\f$ is interpolated on the mesh domain 
using a piecewise linear function which is linear 
in each cell of the Voronoi diagram of the mesh vertices. 

The optimizer `lloyd_optimize_mesh_2()` works in iterative steps. 
At each iteration, mesh vertices are moved into 
positions that bring to zero the energy gradient 
and the Delaunay triangulation is updated. 
Vertices on the mesh boundaries are handled 
in a special way so as to preserve an accurate 
representation of the domain boundaries. 

\tparam CDT is required to be 2D Triangulation.
it provides the initial mesh 
and is modified by the algorithm 
to represent the final optimized mesh. 

The function has an optional `max_iteration_number` parameter,
which defines the number of iterations for the Optimization process.


*/

template <typename CDT> 
  void
  lloyd_optimize_mesh_2(CDT& cdt,
                        int& max_iteration_number=0);

} /* namespace CGAL */
