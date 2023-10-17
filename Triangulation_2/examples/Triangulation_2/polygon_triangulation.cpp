#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
//
#include <CGAL/Gmpq.h>
#include <CGAL/Cartesian.h>
//
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/draw_triangulation_2.h>
#include <CGAL/mark_domain_in_triangulation.h>
#include <CGAL/Polygon_2.h>

#include <iostream>
#include <unordered_map>

#include <boost/property_map/property_map.hpp>

typedef CGAL::Gmpq NT3;
typedef CGAL::Cartesian<NT3> K;
//typedef CGAL::Exact_predicates_inexact_constructions_kernel       K;
typedef CGAL::Triangulation_vertex_base_2<K>                      Vb;
typedef CGAL::Constrained_triangulation_face_base_2<K>            Fb;
typedef CGAL::Triangulation_data_structure_2<Vb,Fb>               TDS;
typedef CGAL::Exact_predicates_tag                                Itag;
typedef CGAL::Constrained_Delaunay_triangulation_2<K, TDS, Itag>  CDT;
typedef CDT::Face_handle                                          Face_handle;
typedef CDT::Point                                                Point;
typedef CGAL::Polygon_2<K>                                        Polygon_2;

int main( )
{
  //construct two non-intersecting nested polygons
  Polygon_2 polygon1;
  polygon1.push_back(Point(0,0));
  polygon1.push_back(Point(2,0));
  polygon1.push_back(Point(2,2));
  polygon1.push_back(Point(1,1.75));
  polygon1.push_back(Point(0,2));
  Polygon_2 polygon2;
  polygon2.push_back(Point(0.5,0.5));
  polygon2.push_back(Point(1.5,0.5));
  polygon2.push_back(Point(1.5,1.5));
  polygon2.push_back(Point(0.5,1.5));
#if 0
  polygon2.push_back(Point(0.5,-0.5));
  polygon2.push_back(Point(1.5,-0.5));
  polygon2.push_back(Point(1.5,0.5));
  polygon2.push_back(Point(0.5,0.5));
#endif
  //Insert the polygons into a constrained triangulation
  CDT cdt;
  cdt.insert_constraint(polygon1.vertices_begin(), polygon1.vertices_end(), true);
  cdt.insert_constraint(polygon2.vertices_begin(), polygon2.vertices_end(), true);
  //cdt.insert_constraint(Point(0.25, 0.25), Point(0.25, 1.75));

  std::unordered_map<Face_handle, bool> in_domain_map;
  boost::associative_property_map< std::unordered_map<Face_handle,bool> >
    in_domain(in_domain_map);

  //Mark facets that are inside the domain bounded by the polygon
  CGAL::mark_domain_in_triangulation(cdt, in_domain);

  std::map<CDT::Vertex_handle, int> vtxmap;
  int cntr = 0;
  for(CDT::Vertex_handle vv : cdt.finite_vertex_handles())
  {
      auto xx = (*vv);
      const CDT::Point &p = cdt.point(vv);
      //std::cout << "vh: " << xx << std::endl;
      std::cout << cntr << " : " << p.x() << " " << p.y() << std::endl;
      vtxmap.insert(std::make_pair(vv, cntr));
      cntr++;
  }
  //
  unsigned int count=0;
  for (Face_handle f : cdt.finite_face_handles())
  {
      if ( get(in_domain, f) ) {
          auto v0 = (*f).vertex(0);
          auto v1 = (*f).vertex(1);
          auto v2 = (*f).vertex(2);
          std::cout << "f " << count << " : ";
          std::cout << vtxmap[v0] << " ";
          std::cout << vtxmap[v1] << " ";
          std::cout << vtxmap[v2] << std::endl;
          ++count;
      }
  }
  //cdt.infinite_vertex()

  std::cout << "There are " << count << " faces in the domain." << std::endl;
  std::cout << "all : " << cdt.number_of_faces() << std::endl;
  assert(count > 0);
  assert(count < cdt.number_of_faces());

  CGAL::draw(cdt, in_domain);
  return 0;
}
