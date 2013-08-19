// Copyright (c) 2009 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
//
//
// Author(s)     : Raul Gallegos, Jane Tournois
//
//******************************************************************************
// File Description : Odt move function
//******************************************************************************

#ifndef CGAL_ODT_OPTIMIZE_MESH_2_H
#define CGAL_ODT_OPTIMIZE_MESH_2_H


#include <CGAL/Mesh_2/Mesh_global_optimizer.h>
#include <CGAL/Mesh_2/Odt_move.h>
#include <CGAL/Mesh_2/Mesh_sizing_field.h>

namespace CGAL
{
  template <typename CDT> 
  void
  odt_optimize_mesh_2(CDT& cdt,
                      int& max_iteration_number)
{
  typedef Mesh_2::Mesh_sizing_field<CDT>                    Sizing;
  typedef typename Mesh_2::Odt_move<CDT,Sizing>             Move;
  typedef typename Mesh_2::Mesh_global_optimizer<CDT,Move>  Odt_optimizer;
  
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