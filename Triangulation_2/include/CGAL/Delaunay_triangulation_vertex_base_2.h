#ifndef CGAL_DELAUNAY_TRIANGULATION_VERTEX_BASE_2_H
#define CGAL_DELAUNAY_TRIANGULATION_VERTEX_BASE_2_H

#include <CGAL/Triangulation_vertex_base_2.h>

namespace CGAL {

template < typename GT,
           typename Vb = Triangulation_vertex_base_2<> >
class Delaunay_triangulation_vertex_base_2 : public Vb
{
public:
	typedef Vb Base;
	typedef typename Gt::Point_2 Point;
	int m_tag;

	template < typename TDS2 >
	struct Rebind_TDS {
		typedef typename Base::template Rebind_TDS<TDS2>::Other Vb2;
		typedef Delaunay_triangulation_vertex_base_2<Gt,Vb2>    Other;
};

private:

public:
	My_vertex_base()
		: Base()
	{
	}
	My_vertex_base(const Point & p, void* f)
		: Base(p,f)
	{
	}
	My_vertex_base(const Point & p)
		: Base(p)
	{
	}

	~My_vertex_base()
	{
		this->set_point(Point(-1000,-1000));
	}

	const int tag() const { return m_tag; }
	int& tag() { return m_tag; }

};
} //namespace CGAL

#endif //CGAL_DELAUNAY_TRIANGULATION_VERTEX_BASE_2_H