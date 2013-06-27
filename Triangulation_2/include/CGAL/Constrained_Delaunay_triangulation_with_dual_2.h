#ifndef CGAL_CONSTRAINED_DELAUNAY_TRIANGULATION_WITH_DUAL_2_H
#define CGAL_CONSTRAINED_DELAUNAY_TRIANGULATION_WITH_DUAL_2_H

#include <CGAL/Constrained_Delaunay_triangulation_2.h>

namespace CGAL {

namespace internal{
template <class Cdt>
class Bvd_cell : public std::pair<typename Cdt::Polygon,			  // cell 
								typename Cdt::Vertex_handle> // generator
{
	typedef typename Bvd_cell<Cdt> Cell;

public:
	Cell()
		: std::pair<typename Cdt::Polygon,	
	 		  			  typename Cdt::Vertex_handle>()
	{
	}
	Cell(typename Cdt::Polygon poly,	
			 typename Cdt::Vertex_handle v)
			 : std::pair<typename Cdt::Polygon,	
			  					 typename Cdt::Vertex_handle>(poly, v)
	{
	}
									
	bool operator<(const Cell& cell) const
	{
		typename Cdt::Point p1 = *((*this).first.left_vertex());
		typename Cdt::Point p2 = *(cell.first.left_vertex());
		return( p1.x() <  p2.x()
			  || (p1.x() == p2.x()  &&  p1.y() < p2.y()));
	};

	typename Cdt::Vertex_handle get_generator() const
	{
		return (*this).second;
	}

	typename Cdt::Polygon get_polygon() const
	{
		return (*this).first;
	}
}; //end CLASS Bvd_cell

} } //namespace CGAL::internal


namespace CGAL {

template <class Gt,
			class Tds = Triangulation_data_structure_2 <
                      Delaunay_triangulation_vertex_base_2<Gt>,
		      Constrained_Delaunay_triangulation_face_base_2<Gt> >,
	  class Itag = No_intersection_tag >
class Constrained_Delaunay_triangulation_with_dual_2 : public Constrained_Delaunay_triangulation_2 //public CGAL::Constrained_triangulation_plus_2<Cdt> 
{
public:
	// typedefs for basic primitives
	
  typedef Constrained_Delaunay_triangulation_2<Gt,Tds,Itag>    Ctr;
  typedef Constrained_Delaunay_triangulation_with_dual_2 Cdt;
  typedef typename Ctr::Geom_traits      Geom_traits;
  typedef typename Ctr::Intersection_tag Intersection_tag;

  typedef typename Ctr::Constraint    Constraint;
  typedef typename Ctr::Vertex_handle Vertex_handle;
  typedef typename Ctr::Face_handle   Face_handle;
  typedef typename Ctr::Edge          Edge;
  typedef typename Ctr::Finite_faces_iterator Finite_faces_iterator;
  typedef typename Ctr::Face_circulator       Face_circulator;
  typedef typename Ctr::size_type             size_type;
  typedef typename Ctr::Locate_type           Locate_type;
 
  typedef typename Ctr::List_edges List_edges;  
  typedef typename Ctr::List_faces List_faces;
  typedef typename Ctr::List_vertices  List_vertices;
  typedef typename Ctr::List_constraints List_constraints;
  typedef typename Ctr::Less_edge less_edge;
  typedef typename Ctr::Edge_set Edge_set;

	typedef typename Cdt::Geom_traits kernel;
	typedef typename kernel::FT FT;
	typedef typename kernel::Point_2 Point;
	typedef typename Cdt::Segment  Segment;
	typedef typename kernel::Vector_2  Vector;
	typedef typename Cdt::Triangle Triangle;
	typedef CGAL::Polygon_2<kernel,std::vector<Point> > Polygon;
	typedef typename Cdt::Edge Edge;

	// typedefs for triangulation
	typedef typename CDT::Edge_iterator   Edge_iterator;
	typedef typename CDT::Vertex_iterator Vertex_iterator;
	typedef typename CDT::Edge_circulator Edge_circulator;

	// typedefs for BVD
	typedef CGAL::Polygon_2<kernel, std::vector<Point> > Polygon;
	typedef Bvd_cell<CDT> Bvd_cell;
	typedef typename std::list<Bvd_cell> Bvd; 

protected:
		FT m_bounding_box[4]; // xmin,xmax,ymin,ymax
		Bvd m_bvd;

private:
		FT m_color_table[300]; //3 numbers for each color (100 colors)

public:
		CCDT() 
		{
			m_bounding_box[0] = m_bounding_box[2] = 0.0;
			m_bounding_box[1] = m_bounding_box[3] = 1.0;
			this->construct_bvd();

			this->contruct_color_table();
		}

		virtual ~CCDT() 
		{
			m_bvd.clear();
		}

	enum {INSIDE = -1,
				UNDETERMINED = 0,
				OUTSIDE = 1};


	//----------------------------------------------------------------
	//--------------------ABOUT FACES SIGHT---------------------------
	//----------------------------------------------------------------

public:

	// blind = false IFF each face sees its circumcenter
	void tag_all_faces_blind(const bool blind) 
	{
		All_faces_iterator f = NULL;
		for(f = this->all_faces_begin();
				f != this->all_faces_end();
				f++)
			f->blind() = blind;
	}

	// blind test for each face
	// if true, set corresponding barrier constraint
	void tag_faces_blind()
	{
		if(dimension() < 2)
			return;

		tag_all_faces_blind(false);

		// for each constrained edge, mark blinded triangles
		for(Finite_edges_iterator e = this->finite_edges_begin();
				e != this->finite_edges_end();
				++e)
		{
			Edge edge = *e;
			if(this->is_constrained(edge))
			{
				tag_neighbors_blind(edge.first, edge);
				Edge twin = this->twin_edge(edge);
				tag_neighbors_blind(twin.first, twin);
			}
		}
	}

	// returns the "other half" of an edge
	// (shared by edge.first and edge.first->neighbor(edge.second))
	Edge twin_edge(Edge e)
	{
		Face_handle neighb = e.first->neighbor(e.second);
		return Edge(neighb,
								neighb->index(e.first));
	}

	// test face for blind with respect to the edge constraint
	void tag_face_blind(Face_handle& f, 
											const Edge constraint)
	{  
		Point c = this->circumcenter(f);

		// les 2 sommets de l'arete (suivant f orientee ds le sens direct)
		Point a = constraint.first->vertex(ccw(constraint.second))->point();  
		Point b = constraint.first->vertex(cw(constraint.second))->point();

		if(segment_hides_circumcenter(this->segment(constraint), this->triangle(f)))	
		{
			f->blind() = true;
			f->blinding_constraint() = constraint;
		}
	}

	// predicate: returns true if the triangle tr and its circumcenter
	// are on the opposite side of the segment seg
	bool segment_hides_circumcenter(const Segment seg,
																	const Triangle tr)
	{
		Point A = seg.source();
		Point B = seg.target();
		double dX = B.x() - A.x();
		double dY = B.y() - A.y();

		Point p0 = tr[0];
		Point p1 = tr[1];
		Point p2 = tr[2];
		double R0 = p0.x()*p0.x() + p0.y()*p0.y();
		double R1 = p1.x()*p1.x() + p1.y()*p1.y();
		double R2 = p2.x()*p2.x() + p2.y()*p2.y();
		double denominator = (p1.x()-p0.x())*(p2.y()-p0.y()) +
			                    (p0.x()-p2.x())*(p1.y()-p0.y());

		double to_test = 2*denominator * (A.x()*dY - A.y()*dX)
										- (R2-R1) * (p0.x()*dX + p0.y()*dY)
										- (R0-R2) * (p1.x()*dX + p1.y()*dY)
										- (R1-R0) * (p2.x()*dX + p2.y()*dY);
		if( to_test > 0 )
			return false;
		else 
			return true;
	}

  // tags with their sights, with respect to the Edge constraint,
	// seed and its neighbor faces, on the same side of Edge than seed.
	void tag_neighbors_blind(Face_handle& seed,
													 const Edge constraint)
	{
		CGAL_assertion(this->is_constrained(constraint));
		if(!this->is_infinite(seed) 
			 && !seed->blind() 
			 && triangle(seed).area() != 0) //to avoid flat triangles outside the domain
		{
			std::stack<Face_handle> faces;
			faces.push(seed);

			while(!faces.empty())
			{
				Face_handle f = faces.top();
				faces.pop();
				this->tag_face_blind(f, constraint);
				if( f->blind())
					this->push_unvisited_neighbors(f, faces);
			}   
		}
	}

	// puts in the stack the unvisited (un-tagged) neighbor faces of f
	void push_unvisited_neighbors(Face_handle& f,
																std::stack<Face_handle>& faces)
	{
		for(int i=0; i<3; i++)
		{
			Face_handle fi = f->neighbor(i);
			Edge edge_i = Edge(f, i);
			if(	!this->is_constrained(edge_i) &&
					!fi->blind() &&
					!is_infinite(fi)) 
				faces.push(fi);
		}
	}   

//--------------------------------------------------------------------
//-----------------------COMPONENTS-----------------------------------
//--------------------------------------------------------------------

	void insert_component(std::list<Point>& component)
	{
		if(component.size() < 3)
			return;
		std::list<Point>::iterator it;
		for(it = component.begin();
			it != component.end();
			it++)
		{
			const Point& p1 = *it;
			Point p2;
			it++;
			if(it == component.end())
				p2 = *component.begin();
			else
				p2 = *it;
			it--;

			insert_constraint(p1,p2);
		}
	}

	void insert_component_no_constraints(std::list<Point>& component)
	{
		if(component.size() < 3)
			return;
		std::list<Point>::iterator it;
		for(it = component.begin();
			it != component.end();
			it++)
		{
			const Point& p = *it;
			insert(p);
		}
	}

	void insert_spay(const Point& p,
		               const unsigned int nb,
									 const double range)
	{
		for(unsigned int i=0;i<nb;i++)
		{
			double rx = (double)rand() / (double)RAND_MAX;
			double ry = (double)rand() / (double)RAND_MAX;
			rx *= rx;
			ry *= ry;
			double dx = (rx - 0.5) * range;
			double dy = (ry - 0.5) * range;

			Point r(p.x() + dx,
				      p.y() + dy);
			insert(r);
		}
	}

	void insert_circle(const Point& center,
										 const double ray,
										 const unsigned int nb_points)
	{
		Vertex_handle v1 = this->insert(Point(center.x()+ray, center.y()));
		Vertex_handle begin = v1;
		Vertex_handle v2;
		for(unsigned int i=1; i<nb_points; i++)
		{
			double angle = (double)i/nb_points * 6.283185307179;
			v2 = this->insert(Point(center.x() + ray * cos(angle),
															center.y() + ray * sin(angle)));
			this->insert_constraint(v1, v2);
			v1 = v2;
		}
		this->insert_constraint(v1,begin);
	}

	// draw inside delaunay edges
	void gl_draw_unconstrained_edges(float line_width,
		unsigned char r,
		unsigned char g,
		unsigned char b)
	{
		::glColor3ub(r,g,b);
		::glLineWidth(line_width);
		::glBegin(GL_LINES);
		Edge_iterator e;
		for(e  = edges_begin(); 
				e != edges_end(); 
				e++) 
		{
			if((*e).first->is_constrained((*e).second))
				continue;

			Face_handle f1 = (*e).first;
			Face_handle f2 = f1->neighbor((*e).second);
			if(f2 == NULL | f1 == NULL)
				continue;

			if(f1->location() == OUTSIDE &&
				 f2->location() == OUTSIDE)
				 continue;

			Point p1 = (*e).first->vertex(ccw((*e).second))->point();
			Point p2 = (*e).first->vertex(cw((*e).second))->point();
			::glVertex2d(CGAL_NTS to_double(p1.x()),CGAL_NTS to_double(p1.y()));
			::glVertex2d(CGAL_NTS to_double(p2.x()),CGAL_NTS to_double(p2.y()));
		}
		::glEnd();
	}

	
	// draw all the delaunay edges (in and out)
	void gl_draw_all_edges(float line_width,
		unsigned char r,
		unsigned char g,
		unsigned char b)
	{
		::glColor3ub(r,g,b);
		::glLineWidth(line_width);
		::glBegin(GL_LINES);
		Edge_iterator e;
		for(e  = edges_begin(); 
				e != edges_end(); 
				e++) 
		{
			Face_handle f1 = (*e).first;
			Face_handle f2 = f1->neighbor((*e).second);
	/*		if(f2 == NULL | f1 == NULL)
				continue;
*/
		/*	if(f1->location() == OUTSIDE &&
				 f2->location() == OUTSIDE)
				 continue;
*/
			Point p1 = (*e).first->vertex(ccw((*e).second))->point();
			Point p2 = (*e).first->vertex(cw((*e).second))->point();
			::glVertex2d(CGAL_NTS to_double(p1.x()),CGAL_NTS to_double(p1.y()));
			::glVertex2d(CGAL_NTS to_double(p2.x()),CGAL_NTS to_double(p2.y()));
		}
		::glEnd();
	}


	Vertex_handle select_vertex(const Point& p)
	{
		Vertex_handle v = nearest_vertex(p, locate(p));
		if(!this->incident_constraints(v))
			return v;
		return NULL;
	}

	Vertex_handle nearest_vertex(const Point& p, 
															 Face_handle& f)
	{
		f = locate(p,f);
	  typename kernel::Compare_distance_2 compare_distance 
			= this->geom_traits().compare_distance_2_object();
		Vertex_handle nn =  !is_infinite(f->vertex(0)) ? f->vertex(0) : f->vertex(1);
	
		if ( !is_infinite(f->vertex(1)) 
			&& compare_distance(p, f->vertex(1)->point(), nn->point()) == CGAL::SMALLER) 
			nn=f->vertex(1);
		if ( !is_infinite(f->vertex(2)) 
			&& compare_distance(p, f->vertex(2)->point(), nn->point()) == CGAL::SMALLER) 
			nn=f->vertex(2);

		look_nearest_neighbor(p,f,0,nn);
		look_nearest_neighbor(p,f,1,nn);
		look_nearest_neighbor(p,f,2,nn);
		return nn;
	}

	void look_nearest_neighbor(const Point& p,
														 Face_handle f,
														 int i,
														 Vertex_handle& nn) const
	{
		Face_handle  ni = f->neighbor(i);
		if ( CGAL::ON_POSITIVE_SIDE != side_of_oriented_circle(ni,p) ) return;

		typename kernel::Compare_distance_2 compare_distance 
			= this->geom_traits().compare_distance_2_object();
		i = ni->index(f);
		if ( !is_infinite(ni->vertex(i)) 
			   & compare_distance(p, ni->vertex(i)->point(), nn->point())  == CGAL::SMALLER)  
			nn=ni->vertex(i);
    
		// recursive exploration of triangles whose circumcircle contains p
		look_nearest_neighbor(p, ni, ccw(i), nn);
		look_nearest_neighbor(p, ni, cw(i), nn);
	} 


	// draw constrained edges
	void gl_draw_constrained_edges(float line_width,
		unsigned char r,
		unsigned char g,
		unsigned char b)
	{
		::glColor3ub(r,g,b);
		::glLineWidth(line_width);
		::glBegin(GL_LINES);
		Edge_iterator e;
		for(e  = edges_begin(); 
	 			e != edges_end(); 
				e++) 
		{
			if((*e).first->is_constrained((*e).second))
			{
				Point p1 = (*e).first->vertex(ccw((*e).second))->point();
				Point p2 = (*e).first->vertex(cw((*e).second))->point();
				::glVertex2d(CGAL_NTS to_double(p1.x()),CGAL_NTS to_double(p1.y()));
				::glVertex2d(CGAL_NTS to_double(p2.x()),CGAL_NTS to_double(p2.y()));
			}
		}
		::glEnd();
	}


	// draw inside face
	void gl_draw_inside_faces(unsigned char r,
		                        unsigned char g,
		                        unsigned char b)
	{
		::glColor3ub(r,g,b);
		::glBegin(GL_TRIANGLES);
		for(Finite_faces_iterator f = finite_faces_begin();
			f != finite_faces_end();
			f++)
		{
			if(f->location() == INSIDE)
			{
				const Point& p1 = f->vertex(0)->point();
				const Point& p2 = f->vertex(1)->point();
				const Point& p3 = f->vertex(2)->point();
				::glVertex2d(CGAL_NTS to_double(p1.x()),CGAL_NTS to_double(p1.y()));
				::glVertex2d(CGAL_NTS to_double(p2.x()),CGAL_NTS to_double(p2.y()));
				::glVertex2d(CGAL_NTS to_double(p3.x()),CGAL_NTS to_double(p3.y()));
			}
		}
		::glEnd();
	}

	void tag_faces_location(const int location)
	{
		for(All_faces_iterator f = all_faces_begin();
			f != all_faces_end();
			f++)
			f->location() = location;
	}


	/*--------------------------------------------------------------
	------------------- BVD CELLS CONSTRUCTION -------------------
	--------------------------------------------------------------*/

	void construct_bvd()
	{
		m_bvd.clear();
		
	    tag_faces_blind();

		for(Finite_vertices_iterator v = this->finite_vertices_begin();
				v != this->finite_vertices_end();
				++v)
		{
			if(!this->cell_is_infinite(v))
				m_bvd.push_back(this->cell(v));
		}
	}

	void update_bvd()
	{
		this->construct_bvd();
	}

	// assemble a cell of the bounded Voronoi diagram
	// incident to vertex v
	typename Bvd_cell cell(Vertex_handle v)
	{
		Polygon polygon;

		ASSERT(!is_infinite(v));
		Face_circulator face = this->incident_faces(v);
		Face_circulator begin = face;
		Face_circulator next = face;

		Point intersection;
		Line line; 

		CGAL_For_all(face, begin)
		{
			next++;
			line = Line(this->circumcenter(face), this->circumcenter(next));

			if(!face->blind()) //face sees
			{
				polygon.push_back(this->circumcenter(face));
				if(next->blind())  //next doesn't
				{
					CGAL_assertion(do_intersect(line, this->segment(next->blinding_constraint())));
					assign(intersection, CGAL::intersection(line, Line(this->segment(next->blinding_constraint()))));
					polygon.push_back(intersection);	  
				}
			}
			else //face doesn't see
			{
				if(!next->blind()) //next sees
				{
					CGAL_assertion(do_intersect(line, this->segment(face->blinding_constraint())));
					assign(intersection, CGAL::intersection(line, Line(this->segment(face->blinding_constraint()))));
					polygon.push_back(intersection);		  
				}
				else //next doesn't
				{
					if(	 face->blinding_constraint() != next->blinding_constraint()
						&& face->blinding_constraint() != this->twin_edge(next->blinding_constraint())) // the 2 blinding_constraints are different
					{
						CGAL_assertion(do_intersect(line, this->segment(face->blinding_constraint())));
						assign(intersection, CGAL::intersection(line, Line(this->segment(face->blinding_constraint()))));
						polygon.push_back(intersection);

						Point intersection2;
						CGAL_assertion(do_intersect(line, this->segment(next->blinding_constraint())));
						assign(intersection2, CGAL::intersection(line, Line(this->segment(next->blinding_constraint()))));
						polygon.push_back(intersection2);
					}
					//else: it's the same constraint--> do nothing
				}
			}
		}//end CGAL_For_all

		return Bvd_cell(polygon, v); 
	}  


	//returns true IFF generator's cell is clipped by at least one constrained Edge
	bool cell_is_clipped(const Vertex_handle generator)
	{
		Face_circulator face = this->incident_faces(generator);
		Face_circulator begin = face;
		CGAL_For_all(face, begin){
			if(face->blind())
				return true;
		}
		return false;
	}

	//returns true IFF generators's cell is finite
	bool cell_is_infinite(const Vertex_handle generator)
	{
		Face_circulator face = this->incident_faces(generator);
		Face_circulator begin = face;
		CGAL_For_all(face, begin){
			if(this->is_infinite(face))
				return true;
		}
		return false;
	}

	

//----------------------------------------------------
//------------------ TOOLS ---------------------------
//----------------------------------------------------

	//removes all the unconstrained vertices of the mesh
	void clear_unconstrained_vertices()
	{
		for(Finite_vertices_iterator fvt = finite_vertices_begin();
				fvt != finite_vertices_end();
				fvt++)
		{
			if(!this->are_there_incident_constraints(fvt))
				this->remove(fvt);
		}
	}


	void cdt_remove_constraint(Vertex_handle va, Vertex_handle vb)
	{
		Face_handle f;  
    int index;
		bool b = this->is_edge(va, vb, f, index);
		VERIFY(b);  
		this->remove_constraint(f,index);
	}


	bool is_inside(const Point& query)
	{
		Face_handle f = locate(query);
		if(f == NULL)
			return false;
		if(f->location() == INSIDE)
			return true;
		return false;
	}

	bool incident_constraints(Vertex_handle v)
	{
		Edge_circulator e = incident_edges(v);
		Edge_circulator end = e;
		CGAL_For_all(e,end)
			if((*e).first->is_constrained((*e).second))
				return true;
		return false;
	}


	//---------------------------------------------------------------
	//------------------- DRAW BVD ----------------------------------
	//---------------------------------------------------------------

	// draw Voronoi edges
	// note: Voronoi of a constrained Delaunay
	// triangulation does not have any precise definition.
	void gl_draw_bounded_voronoi_diagram()
	{
		if(dimension() < 2)
			return;

		::glBegin(GL_LINES);

		for(Edge_iterator e = edges_begin();
  			e != edges_end();
	  		++e)
		{
			Edge edge = *e;

			// les 2 faces adjacentes a edge
			Face_handle f1 = edge.first;
			Face_handle f2 = f1->neighbor(edge.second); 

			if(f1 == NULL || f2 == NULL)
				continue;

			// 2 finite faces
			if( !is_infinite(f1) && !is_infinite(f2))
				if( !f1->blind() && !f2->blind())
					gl_draw_voronoi_edge(f1, f2, edge);
				else
					gl_draw_clipped_voronoi_edge(f1, f2, edge);

			//une seule face infinie, l'autre finie
			else if( !is_infinite(f1) && is_infinite(f2) && triangle(f1).area() != 0)
				gl_draw_ray_bounded_voronoi_diagram(edge, f1);
			else if( is_infinite(f1) && !is_infinite(f2) && triangle(f2).area() != 0)
				gl_draw_ray_bounded_voronoi_diagram(edge, f2);

			// les 2 faces infinies <=> pslg de dimension 1
		}  
		::glEnd();
	}

	//f1 and f2 can't be of area 0 simultaneously
	//(see the use of the function above)
	void gl_draw_voronoi_edge(const Face_handle f1,
														const Face_handle f2,
														const Edge& e)
	{
		if(triangle(f1).area() == 0 && triangle(f2).area() == 0)
			return;
		else if(triangle(f2).area()!=0 && triangle(f1).area()==0)
			gl_draw_ray_bounded_voronoi_diagram(e,f2);
		else if(triangle(f1).area()!=0 && triangle(f2).area()==0)
			gl_draw_ray_bounded_voronoi_diagram(e,f1);
		else
		{
			Point a = circumcenter(f1);
			Point b = circumcenter(f2);
			::glVertex2d(a.x(),a.y());
			::glVertex2d(b.x(),b.y());
		}
	}

	void debug()
	{
		K::Do_intersect_2 do_intersect = K().do_intersect_2_object();
		K::Intersect_2 intersect = K().intersect_2_object();

		// UN CAS QUI NE MARCHE PAS
		Point c3 = Point( -0.595731, 0.659805);
		Point c4 = Point(-0.780416, 0.563209);
		Line l2 = Line(c3, c4);

		Segment s2 = Segment(Point(-0.686405, 0.858313),
			Point(-0.686405, -0.45309));

		CGAL::Object obj = intersect(l2, s2);
		if(obj.is_empty())
			int bidon = 0;
		//std::type_info type_obj = obj.type().name();

		Point p;
		Segment s;
		if (assign(p, obj)) 
			std::cout<< "2nd point OK: " <<p<<std::endl;
		else
			if(assign(s, obj))
			{
				int bidon = 0;
			}
			else
			{
				if(do_intersect(l2, s2))
				{
					int bidon = 0;
				}
			}
	}

	void gl_draw_debug()
	{
		Point a = Point( -0.595731, 0.659805);
		Point b = Point(-0.780416, 0.563209);
		Point c = Point( -0.686405, 0.858313);
		Point d = Point(-0.686405, -0.45309);
		glBegin(GL_LINES);
		glVertex2d(a.x(),a.y());
		glVertex2d(b.x(),b.y());
		glVertex2d(c.x(),c.y());
		glVertex2d(d.x(),d.y());
		glEnd();
	}

	void gl_draw_clipped_voronoi_edge(const Face_handle f1,
																		const Face_handle f2, 
																		const Edge& e)
	{
		Point c1 = circumcenter(f1);
		Point c2 = circumcenter(f2);
		Point intersection;
		Line c1c2 = Line(c1,c2);
		
		// f1 voit, pas f2  (f1 peut etre d'aire 0)
		if(!f1->blind() && f2->blind())  
		{
			if(triangle(f1).area()!=0)
			{
				VERIFY(assign(intersection, CGAL::intersection(c1c2, Line(segment(f2->blinding_constraint())))));
				::glVertex2d(c1.x(),c1.y());
				::glVertex2d(intersection.x(),intersection.y());
			}else
				gl_draw_ray_bounded_voronoi_diagram(e,f1);
		}
		// f2 voit, pas f1  (f2 peut etre d'aire 0)
		else if(!f2->blind() && f1->blind())  
		{
			if(triangle(f2).area()!=0)
			{
			  VERIFY(assign(intersection, CGAL::intersection(c1c2, Line(segment(f1->blinding_constraint())))));
				::glVertex2d(c2.x(),c2.y());
				::glVertex2d(intersection.x(),intersection.y());
			}else
				gl_draw_ray_bounded_voronoi_diagram(e,f2);
		}
		// f1 et f2 ne voient pas (aucune ne peut etre d'aire 0)
		else if(f1->blind() && f2->blind()) 
		{
			if( f1->blinding_constraint() != twin_edge(f2->blinding_constraint()) )
			{
				Point intersection2;
				Line l1 = Line(segment(f1->blinding_constraint()));
				Line l2 = Line(segment(f2->blinding_constraint()));

				VERIFY(assign(intersection, CGAL::intersection(c1c2, l1)));
				VERIFY(assign(intersection2, CGAL::intersection(c1c2, l2)));

				::glVertex2d(intersection.x(),intersection.y());
				::glVertex2d(intersection2.x(),intersection2.y());
			}
		}
	}

	void gl_draw_ray_bounded_voronoi_diagram(const Edge& edge,
																					 const Face_handle finite_face)
	{
		Segment s = segment(edge);
		Point c = circumcenter(finite_face);

		Vector e = (edge.first == finite_face) ? (s.target() - s.source()) : (s.source() - s.target());
		Vector n = Vector(e.y(), -e.x()); //rotation de -pi/2 => normale sortante de l'enveloppe convexe
		Point target = c + 100.0 * n;

		if(!finite_face->blind())
		{
			::glVertex2d(c.x(),c.y());
			::glVertex2d(target.x(),target.y());
		}
		else
			draw_clipped_voronoi_ray(edge, finite_face, target);
	}

	//draw the ray "clipped" on clipping_edge
	void draw_clipped_voronoi_ray(const Edge& clipping_edge,
																const Face_handle finite_face,
																const Point target)
	{
		Point on_edge;
		Segment clipping_constraint;
		if(finite_face->location() == INSIDE)
		{
			clipping_constraint = segment(clipping_edge);
			on_edge = CGAL::midpoint(clipping_constraint.source(),clipping_constraint.target());	
		}
		else if(finite_face->location() == OUTSIDE)
		{
			clipping_constraint = segment(finite_face->blinding_constraint());
			Line l1 = Line(clipping_constraint.source(), clipping_constraint.target());
			Line l2 = Line(target, circumcenter(finite_face));
			VERIFY(assign(on_edge, CGAL::intersection(l1, l2)));
		}

		::glVertex2d(on_edge.x(),on_edge.y());
		::glVertex2d(target.x(),target.y());
	}

	//typedef typename Cdt::Locate_type Locate_type;
	bool inside_convex_hull(const Point& query)
	{
		// locate
		Locate_type locate_type;
		int li;
		Face_handle f = locate(query,locate_type,li);

		// check inside/outside convex hull of the triangulation
		return (locate_type != this->OUTSIDE_CONVEX_HULL);
	}


	// fills the "blind" faces
	void gl_show_blind_faces()
	{
		this->tag_faces_blind();
		::glBegin(GL_TRIANGLES);
		for(Finite_faces_iterator f = finite_faces_begin();
			f != finite_faces_end();
			++f)
		{
			if(f->blind())
			{
				const Point& a = f->vertex(0)->point();
				const Point& b = f->vertex(1)->point();
				const Point& c = f->vertex(2)->point();
				::glVertex2d(a.x(),a.y());
				::glVertex2d(b.x(),b.y());
				::glVertex2d(c.x(),c.y());
			}
		}     
		::glEnd();
	}


	void gl_draw_constructed_bvd()
	{
		this->construct_bvd();
		int rd_r, rd_g, rd_b; //random color components

		std::list<Bvd_cell>::iterator cell_it;
		unsigned int n = 0;
		for(cell_it = m_bvd.begin();
				cell_it != m_bvd.end();
				cell_it++)
		{
			Polygon poly = cell_it->first;
			::glBegin(GL_TRIANGLES);
			rd_r = (int)m_color_table[n%300]; // couleurs pastel
			rd_g = (int)m_color_table[(n+1)%300];
			rd_b = (int)m_color_table[(n+2)%300];
			::glColor3ub(rd_r,rd_g,rd_b);
			for(unsigned int i = 0; i < poly.size(); i++)
			{
				::glVertex2d(poly.vertex(i).x(), poly.vertex(i).y());
				::glVertex2d(poly.vertex((i+1)%poly.size()).x(), poly.vertex((i+1)%poly.size()).y());
				::glVertex2d(cell_it->second->point().x(), cell_it->second->point().y());
			}
			::glEnd();
			n = n+3;
		}
	}

	int my_random(int min, int max){
		return (int) (min + rand() % (max - min) );
	}

	void contruct_color_table()
	{
		for(unsigned int i=0; i < 300; i++)
			m_color_table[i] = my_random(200,255);
	} 
	

//---------------------------------------------------------------
//--------------------- SAVE MESH -------------------------------
//---------------------------------------------------------------


	//save the cdt and/or the bvd as an eps file
	bool save_eps(char *pFilename, 
								const bool save_cdt,
								const bool save_bvd)
	{
		FILE *pFile = fopen(pFilename,"wt");
		if(pFile == NULL)
			return false;

		fprintf(pFile,"%%!PS-Adobe-2.0 EPSF-2.0\n");
		fprintf(pFile,"%%%%BoundingBox: 0 0 500 500\n");
		fprintf(pFile,"%%%%EndComments\n");
		fprintf(pFile,"gsave\n");
		fprintf(pFile,"0.1 setlinewidth\n");

		// color and thickness macros
		fprintf(pFile,"/color_edges { 0 0 0 setrgbcolor } bind def\n");
		fprintf(pFile,"/linewidth_edges { 0.1 setlinewidth } bind def\n");

		fprintf(pFile,"/color_constraints { 1 0 0 setrgbcolor } bind def\n");
		fprintf(pFile,"/linewidth_constraints { 1.0 setlinewidth } bind def\n");

		fprintf(pFile,"/color_bvd { 0 0 1 setrgbcolor } bind def\n");
		fprintf(pFile,"/linewidth_bvd { 0.9 setlinewidth } bind def\n");

		// edge macro -> E
		fprintf(pFile,"\n%% stroke - x1 y1 x2 y2 E\n");
		fprintf(pFile,"/E {newpath moveto lineto stroke} bind def\n");

		//factor over all coordinates
		const double r = 500;

		// output edges
 	  if(save_cdt)
			save_cdt_eps(pFile, r);

		// input constraints (always)
		save_input_constraints_eps(pFile, r);

		//bvd edges
		if(save_bvd)
			save_bvd_eps(pFile, r);

		fputs("grestore\n\n",pFile);
		fputs("showpage\n",pFile);
		fclose(pFile);
		return true;
	}

	bool save_cdt_eps(FILE *pFile, const double r)
	{
		fprintf(pFile,"\n");
		fprintf(pFile,"color_edges\n");
		fprintf(pFile,"linewidth_edges\n");
		for(Edge_iterator e  = this->edges_begin(); 
				e != this->edges_end(); 
				e++) 
		{
			if(is_constrained(*e))
				continue;

			Face_handle f1 = (*e).first;
			Face_handle f2 = f1->neighbor((*e).second);
			if(f1->location() == OUTSIDE && f2->location() == OUTSIDE)
				continue;
			Point a = (*e).first->vertex(ccw((*e).second))->point();
			Point b = (*e).first->vertex( cw((*e).second))->point();
			fprintf(pFile,"%g %g %g %g E\n",r*a.x(),r*a.y(),r*b.x(),r*b.y());
		}
		return true;
	}

	bool save_input_constraints_eps(FILE *pFile, const double r)
	{
		fprintf(pFile,"\n");
		fprintf(pFile,"color_constraints\n");
		fprintf(pFile,"linewidth_constraints\n");
		for(Edge_iterator e  = this->edges_begin(); 
				e != this->edges_end(); 
				e++) 
		{
			if(!is_constrained(*e))
				continue;

			Face_handle f1 = (*e).first;
			Face_handle f2 = f1->neighbor((*e).second);
			if(f1->location() == OUTSIDE && f2->location() == OUTSIDE)
				continue;

			Point a = (*e).first->vertex(ccw((*e).second))->point();
			Point b = (*e).first->vertex( cw((*e).second))->point();
			fprintf(pFile,"%g %g %g %g E\n",r*a.x(),r*a.y(),r*b.x(),r*b.y());
		}
		return true;
	}

	bool save_bvd_eps(FILE *pFile, const double r)
	{
		fprintf(pFile,"\n");
		fprintf(pFile,"color_bvd\n");
		fprintf(pFile,"linewidth_bvd\n");
		
		//to close the cells of the constrained vertices
		for(Edge_iterator eit = this->edges_begin();
				eit != this->edges_end();
				eit++)
		{      
			Face_handle f1 = eit->first;
			Face_handle f2 = f1->neighbor(eit->second);
			if(f1 == NULL || f2 == NULL)
				continue;

			//2 finite faces -> finite edge(s)
			if(!is_infinite(f1) && !is_infinite(f2))
				if(!f1->blind() && !f2->blind())
					save_voronoi_edge_eps(pFile,f1,f2,*eit,r);
				else
					save_clipped_voronoi_edge_eps(pFile,f1,f2,*eit,r);

			//1 finite face & 1 infinite face
			else if(!is_infinite(f1) && is_infinite(f2))
				save_ray_bvd_eps(pFile, *eit, f1, r);
			else if(is_infinite(f1) && !is_infinite(f2))
				save_ray_bvd_eps(pFile, *eit, f2, r);
		}
		return true;
	}

	bool save_voronoi_edge_eps(FILE *pFile, const Face_handle f1, const Face_handle f2, const Edge& e, const double r)
	{
		if(triangle(f1).area()==0 && triangle(f2).area()==0)
			return false;
		else if(triangle(f1).area()==0 && triangle(f2).area()!=0)
			save_ray_bvd_eps(pFile, e, f2, r);
		else if(triangle(f1).area()!=0 && triangle(f2).area()==0)
			save_ray_bvd_eps(pFile, e, f1, r);
		else
		{
			Point c1 = this->circumcenter(f1);
			Point c2 = this->circumcenter(f2);
			fprintf(pFile,"%g %g %g %g E\n",r*c1.x(),r*c1.y(),r*c2.x(),r*c2.y());
		}
		return true;
	}

	bool save_clipped_voronoi_edge_eps(FILE *pFile, const Face_handle f1, const Face_handle f2, const Edge& e, const double r)
	{
		Point c1 = this->circumcenter(f1);
		Point c2 = this->circumcenter(f2);
		Point intersection;
		Line c1c2 = Line(c1,c2);

		//f1 sees, f2 is blind (area(f1) can be 0)
		if(!f1->blind() && f2->blind())
		{
			if(triangle(f1).area()!=0)
			{
				VERIFY(assign(intersection, CGAL::intersection(c1c2,Line(segment(f2->blinding_constraint())))));
				fprintf(pFile,"%g %g %g %g E\n",r*c1.x(),r*c1.y(),r*intersection.x(),r*intersection.y());
			}else
				save_ray_bvd_eps(pFile, e, f1, r);
		}
		//f2 sees, f1 is blind (area(f2) can be 0)
		else if(f1->blind() && !f2->blind())
		{
			if(triangle(f2).area()!=0)
			{
				VERIFY(assign(intersection, CGAL::intersection(c1c2,Line(segment(f1->blinding_constraint())))));
				fprintf(pFile,"%g %g %g %g E\n",r*c2.x(),r*c2.y(),r*intersection.x(),r*intersection.y());
			}else
				save_ray_bvd_eps(pFile, e, f2, r);
		}
		//f1 blind and f2 blind (area(f1)!=0 && area(f2)!=0)
		else if(f1->blind() && f2->blind())
		{
			if(f1->blinding_constraint() != twin_edge(f2->blinding_constraint()))
			{
				VERIFY(assign(intersection, CGAL::intersection(c1c2,Line(segment(f1->blinding_constraint())))));
				Point intersection2;
				VERIFY(assign(intersection2, CGAL::intersection(c1c2,Line(segment(f2->blinding_constraint())))));
				fprintf(pFile,"%g %g %g %g E\n",r*intersection.x(),r*intersection.y(),r*intersection2.x(),r*intersection2.y());
			}
		}
		return true;
	}

	bool save_ray_bvd_eps(FILE *pFile, const Edge& edge, const Face_handle finite_face, const double r)
	{
		//if(triangle(finite_face).area() > min_interior_triangle_area())

		Segment s = segment(edge);
		Point c = this->circumcenter(finite_face);
		Vector e = (edge.first == finite_face) ? (s.target() - s.source()) : (s.source() - s.target());
		Vector n = Vector(e.y(), -e.x());
		Point m = CGAL::midpoint(segment(edge).source(), segment(edge).target());
		Point target = m + 100.0*n;
			
		if(finite_face->blind()) //finite_face blind => clip the ray on the constraint
		{
			Point t_on_border = bounded_target(m, target);
			fprintf(pFile,"%g %g %g %g E\n",r*m.x(),r*m.y(),r*t_on_border.x(),r*t_on_border.y());
		}
		else // finite_face sees its circumcenter
		{
			Point target = c + 100.0*n;
			Point t_on_border = bounded_target(target, c);
			if(triangle(finite_face).area() > 0.01)
				fprintf(pFile,"%g %g %g %g E\n",r*c.x(),r*c.y(),r*t_on_border.x(),r*t_on_border.y());
			else //flat triangle
				fprintf(pFile,"%g %g %g %g E\n",r*m.x(),r*m.y(),r*t_on_border.x(),r*t_on_border.y());
		}
		return true;
	}

	//to "cut" the ray on the border of the bounding box [0,1]x[0,1]
	Point bounded_target(const Point source,
											 const Point target)
	{
		Point b_target = target;
		Line line = Line(source,target);
		Segment seg = Segment(source,target);

		if(do_intersect(seg, Segment(Point(0.,0.),Point(0.,1.))))
			VERIFY(assign(b_target, CGAL::intersection(line, Line(Point(0.,0.),Point(0.,1.)))));
		else if(do_intersect(seg, Segment(Point(0.,1.),Point(1.,1.))))
			VERIFY(assign(b_target, CGAL::intersection(line, Line(Point(0.,1.),Point(1.,1.)))));
		else if(do_intersect(seg, Segment(Point(1.,1.),Point(1.,0.))))
			VERIFY(assign(b_target, CGAL::intersection(line, Line(Point(1.,1.),Point(1.,0.)))));
		else if(do_intersect(seg, Segment(Point(1.,0.),Point(0.,0.))))
			VERIFY(assign(b_target, CGAL::intersection(line, Line(Point(1.,0.),Point(0.,0.)))));

				/*
		if(!fixed_intersection(b_target,seg,Segment(Point(0.,0.),Point(0.,1.))))
			if(!fixed_intersection(b_target,seg,Segment(Point(0.,1.),Point(1.,1.))))
				if(!fixed_intersection(b_target,seg,Segment(Point(1.,1.),Point(1.,0.))))
					if(!fixed_intersection(b_target,seg,Segment(Point(1.,0.),Point(0.,0.))))
						t_on_border = target;
						*/

		return b_target;
	}

	bool fixed_intersection(Point& p,
													const Segment s1,
													const Segment s2)
	{
		if(!CGAL::do_intersect(s1,s2))
			return false;
		else
		{
			VERIFY(assign(p,CGAL::intersection(Line(s1),Line(s2))));
			return true;
		}
	}


	// save a PSLG as mypslg.edg
	bool save_edg(char *pFilename)
	{
		FILE *pFile = fopen(pFilename,"wt");
		if(pFile == NULL)
			return false;

		unsigned int nb = 0;
		for(Edge_iterator e  = this->edges_begin(); 
				e != this->edges_end(); 
				e++) 
		{
			if(is_constrained(*e))
				nb++;
		}
		// print #constraints
		fprintf(pFile,"%d\n",nb);

		for(Edge_iterator e  = this->edges_begin(); 
			e != this->edges_end(); 
			e++) 
		{
			if(is_constrained(*e))
			{
				Point a = (*e).first->vertex(ccw((*e).second))->point();
				Point b = (*e).first->vertex( cw((*e).second))->point();
				fprintf(pFile,"%g %g %g %g\n",a.x(),a.y(),b.x(),b.y());
			}
		}

		fclose(pFile);
		return true;
	}

	bool save_obj(char *pFilename)
	{
		FILE *f = fopen(pFilename,"wt");
		if(f == NULL)
			return false;

		// output vertices
		for(Finite_vertices_iterator v = finite_vertices_begin();
			v != finite_vertices_end();
			v++)
		{
			const double x = CGAL_NTS to_double(v->point().x());
			const double y = CGAL_NTS to_double(v->point().y());
			fprintf(f,"v %f %f 0\n",x,y);
		}

		// output faces
		tag_vertices();
		for(Finite_faces_iterator it = finite_faces_begin();
			it != finite_faces_end();
			it++)
		{
			if(it->location() == INSIDE)
			{
				// 1-based
				const int& i1 = it->vertex(0)->tag() + 1;
				const int& i2 = it->vertex(1)->tag() + 1;
				const int& i3 = it->vertex(2)->tag() + 1;
				fprintf(f,"f %d %d %d\n",i1,i2,i3);
			}
		}

		fclose(f);
		return true;
	}


	void tag_vertices()
	{
		int index = 0;
		for(Finite_vertices_iterator v = finite_vertices_begin();
			v != finite_vertices_end();
			v++)
			v->tag() = index++;
	}




	bool farthest_voronoi_point_from(Vertex_handle v,
																	 Point& pole,
																	 Face_handle& face_pole)
	{
		const Point& p1 = v->point();

		Face_circulator f = incident_faces(v);
		Face_circulator begin = f;
		FT max_sq_distance = 0;
		bool success = false;
		CGAL_For_all(f,begin)
		{
			if(is_infinite(f))
				continue;
			Point p2 = f->circumcenter();
			FT sq_distance = CGAL::squared_distance(p1,p2);
			if(f == begin || 
				(sq_distance > max_sq_distance))
			{
				pole = p2;
				face_pole = f;
				max_sq_distance = sq_distance;
				success = true;
			}
		}
		return success;
	}

	void tag_faces_inside_outside(std::list<Point>& seeds)
	{
		// reset all tags undetermined
		tag_faces_location(INSIDE);

		// pick infinite seed face
		Face_handle seed = infinite_vertex()->face();
		seed->location() = OUTSIDE;
		std::stack<Face_handle> faces;
		faces.push(seed);

		// locate from all seeds
		std::list<Point>::iterator it;
		for(it = seeds.begin();
			it != seeds.end();
			it++)
		{
			const Point& p = *it;
			Face_handle f = locate(p);
			if(f != NULL)
				faces.push(f);
		}

		while(!faces.empty())
		{
			Face_handle f = faces.top();
			faces.pop();
			const int& location = f->location();
			for(unsigned int i=0;i<3;i++)
			{
				if(f->neighbor(i) != NULL)
					if(f->neighbor(i)->location() == INSIDE && 
						!f->is_constrained(i))
					{
						f->neighbor(i)->location() = OUTSIDE;
						faces.push(f->neighbor(i));
					}
			}
		}
	}


	void tag_faces_inside_outside()
	{
		// reset all tags undetermined
		tag_faces_location(UNDETERMINED);

		// pick one seed face
		Face_handle seed = infinite_vertex()->face();
		seed->location() = OUTSIDE;
		std::stack<Face_handle> faces;
		faces.push(seed);

		while(!faces.empty())
		{
			Face_handle f = faces.top();
			faces.pop();
			const int& location = f->location();
			ASSERT(location == OUTSIDE || location == INSIDE);
			for(unsigned int i=0;i<3;i++)
			{
				if(f->neighbor(i) != NULL)
					if(f->neighbor(i)->location() == UNDETERMINED)
					{
						// swap inside / outside if crosses a tagged facet
						bool constrained = f->is_constrained(i);
						f->neighbor(i)->location() = constrained ? -location : location;
						faces.push(f->neighbor(i));
					}
			}
		}
	}

	bool inside(const Point& query)
	{
		// locate
		Locate_type locate_type;
		int li;
		Face_handle f = locate(query,locate_type,li);

		// check inside/outside
		if(locate_type != OUTSIDE_CONVEX_HULL)
			if(f->location() == INSIDE)
				return true;
		return false;
	}

	// draw generators as dots
	void gl_draw_generators(float point_size,
		unsigned char r,
		unsigned char g,
		unsigned char b)
	{
		::glColor3ub(r,g,b);
		::glPointSize(point_size);
		::glBegin(GL_POINTS);
		Point_iterator it;
		for(it  = points_begin(); 
			it != points_end(); 
			it++) 
		{
			const Point& p = *it;
			::glVertex2d(CGAL_NTS to_double(p.x()),
				CGAL_NTS to_double(p.y()));
		}
		::glEnd();
	}

	// draw generators as dots
	void gl_draw_generators_as_discs(unsigned char r,
		unsigned char g,
		unsigned char b,
		const FT ratio = 0.25)
	{
		::glColor3ub(r,g,b);
		GLUquadricObj* pQuadric = ::gluNewQuadric();
		Finite_vertices_iterator v;
		for(v  = finite_vertices_begin(); 
			v != finite_vertices_end(); 
			v++) 
		{
			FT radius = ratio * min_len_incident_edge(v);
			if(radius == 0) // one unique generator
				continue;
			const Point& p = v->point();
			::glPushMatrix();
			::glTranslated(p.x(),p.y(),0);
			::gluDisk(pQuadric,0,radius,60,1); 
			::glPopMatrix();
		}
		::gluDeleteQuadric(pQuadric);
	}

	FT min_len_incident_edge(Vertex_handle v)
	{
		FT min_len = 0;
		Edge_circulator edge = incident_edges(v);
		Edge_circulator begin = edge;
		CGAL_For_all(edge,begin)
		{
			// skip infinite edges
			if(is_infinite(edge))
				continue;
			Point a = (*edge).first->vertex(ccw((*edge).second))->point();
			Point b = (*edge).first->vertex(cw((*edge).second))->point();
			FT d = distance(a,b);
			if(edge == begin || d < min_len)
				min_len = d;
		}
		return min_len;
	}

	FT distance(const Point& a,
		const Point& b)
	{
		return (FT)std::sqrt((a-b)*(a-b));
	}


	// read PSLG (edg format)
  bool read_pslg(const char *pFilename)
  {
		TRACE("read .edg file %s\n",pFilename);
    FILE *pFile = fopen(pFilename,"rt");
    if(pFile == NULL)
      return false;

		unsigned int nb_edges = 0;
		if(fscanf(pFile,"%d",&nb_edges) != 1)
		{
	    fclose(pFile);
			return false;
		}
		TRACE("%d edges\n",nb_edges);

		for(unsigned int i=0;i<nb_edges;i++)
		{
			float x1,y1,x2,y2;
			if(fscanf(pFile,"%f %f %f %f",&x1,&y1,&x2,&y2) != 4)
			{
		    fclose(pFile);
				return false;
			}
			Point a(x1,y1);
			Point b(x2,y2);
			insert_constraint(a,b);
		}
    fclose(pFile);
    return true;
  }


	// compute bounding box
  void compute_bounding_box()
	 {
			Finite_vertices_iterator v;
			for(v = this->finite_vertices_begin();
				  v != this->finite_vertices_end();
					v++)
			{
				if(v == this->finite_vertices_begin())
				{
					m_bounding_box[0] = m_bounding_box[1] = v->point().x();
					m_bounding_box[2] = m_bounding_box[3] = v->point().y();
				}
				else
				{
					m_bounding_box[0] = std::min(m_bounding_box[0],v->point().x());
					m_bounding_box[1] = std::max(m_bounding_box[1],v->point().x());
					m_bounding_box[2] = std::min(m_bounding_box[2],v->point().y());
					m_bounding_box[3] = std::max(m_bounding_box[3],v->point().y());
				}
			}
	 }
	 bool valid_bounding_box()
	 {
		 return (m_bounding_box[1] - m_bounding_box[0]) > 0.0 &&
			      (m_bounding_box[3] - m_bounding_box[2]) > 0.0;
	 }



};// class CCDT

} //namespace CGAL
#endif // CGAL_CONSTRAINED_DELAUNAY_TRIANGULATION_WITH_DUAL_2_H