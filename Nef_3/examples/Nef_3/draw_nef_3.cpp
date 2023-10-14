#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
//
#include <CGAL/Gmpq.h>
#include <CGAL/Cartesian.h>
//
#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/draw_nef_3.h>

#include <CGAL/OFF_to_nef_3.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/convert_nef_polyhedron_to_polygon_mesh.h>

#include <fstream>
#include <iostream>

// double //
//using NT3 = double;
//using CGAL_Kernel3 = CGAL::Cartesian<NT3>;
// Gmpq //
using NT3 = CGAL::Gmpq;
using CGAL_Kernel3 = CGAL::Cartesian<NT3>;
// ???? //
//using NT3 = CGAL::Exact_predicates_exact_constructions_kernel::FT;
//using CGAL_Kernel3 = CGAL::Exact_predicates_exact_constructions_kernel;
using CGAL_Nef_polyhedron3 = CGAL::Nef_polyhedron_3<CGAL_Kernel3>;
using CGAL_Aff_transformation = CGAL_Nef_polyhedron3::Aff_transformation_3;
using CGAL_Polyhedron = CGAL::Polyhedron_3<CGAL_Kernel3>;

//typedef CGAL::Exact_predicates_exact_constructions_kernel   Kernel;
//typedef CGAL::Polyhedron_3<Kernel>                          Polyhedron;
//typedef CGAL::Nef_polyhedron_3<Kernel>                      Nef_polyhedron;

typedef CGAL_Kernel3 Kernel;
typedef CGAL_Polyhedron Polyhedron;
typedef CGAL_Nef_polyhedron3 Nef_polyhedron;

typedef CGAL::Surface_mesh<CGAL::Cartesian<double>::Point_3> Surface_mesh_;
//typedef CGAL::Surface_mesh<Kernel::Point_3> Surface_mesh_;

using namespace CGAL;

template <class Nef_3>
std::size_t
HOGE (std::istream &i_st, Nef_3 &nef_union, bool verb=true)
{
   // Nef
   typedef typename Nef_3::Kernel          Kernel;
   typedef typename Nef_3::Point_3         Point_3;
   typedef typename std::vector<Point_3>   Point_set;
   CGAL::Nef_nary_union_3<Nef_3>            nary_union;
   // input data structure
   typedef double                           Scan_NT;
   typedef CGAL::Cartesian<Scan_NT>        Scan_kernel;
   typedef Scan_kernel::Point_3            Scan_point;
   typedef std::vector<Scan_point>         Scan_point_set;
   typedef Scan_kernel::Vector_3           Scan_vector;

   typedef CGAL::Scanner_OFF<Scan_kernel>  Scan_OFF;
   typedef Scan_OFF::Vertex_iterator       Scan_vertex_it;
   typedef Scan_OFF::Facet_iterator        Scan_facet_it;
   typedef Scan_OFF::Index_iterator        Scan_index_it;

   typedef typename Kernel::Kernel_tag Kernel_tag;
   typedef typename CGAL::number_type_converter_nef_3<Kernel_tag,Kernel> ntc;

   // declarations and defaults
   std::size_t discarded_facets=0;
   std::size_t idx;

   // input: description of polyhedron in object file format (OFF)
   // with cartesian double coordinates
   Scan_OFF scan (i_st);
   std::size_t NOV = scan.size_of_vertices();

   // read and store vertices
   Scan_point_set V_scan;
   V_scan.reserve (NOV);
   Point_set V;
   V.reserve (NOV);

   {
       Scan_point sp;
       std::cout << "x : " << sp.x() << std::endl;
       std::cout << "y : " << sp[1] << std::endl;
   }
   {
       Scan_point sp(1.0,2.0,3.0);
       std::cout << "x : " << sp.x() << std::endl;
       std::cout << "y : " << sp[1] << std::endl;
   }

   Scan_vertex_it v_it = scan.vertices_begin();
   for (idx=0; v_it != scan.vertices_end(); ++v_it, ++idx)
   {  V_scan.push_back (*v_it);
     V.push_back (ntc::convert(*v_it));

     Scan_point &sp = V_scan.back();
     std::cout << sp << std::endl;
     Point_3 &p = V.back();
     std::cout << p << std::endl;
   }
   CGAL_warning ( idx==NOV );
   NOV = idx;

   // for each facet
   Scan_facet_it f_it = scan.facets_begin();
   for (idx=0; f_it != scan.facets_end(); ++f_it, ++idx)
   {  // read facet
      Scan_facet_it::indices_size_type NOI=f_it.size_of_indices(), jdx;
      Scan_point_set V_f_scan;
      V_f_scan.reserve(NOI);
      Point_set V_f;
      V_f.reserve(NOI);

      Scan_index_it ind_it = f_it->begin();
      for (jdx=0; ind_it != f_it->end(); ++ind_it, ++jdx)
      {  // assertion: index out of range?
         CGAL_assertion (*ind_it < NOV );
         V_f_scan.push_back (V_scan[*ind_it]);
         V_f.push_back (V[*ind_it]);
      }
      CGAL_warning ( jdx==NOI );
      NOI = jdx;

      bool is_nef = false;
      CGAL_assertion_msg( V_f.size() >= 1 || !verb, "empty vertex cycle");
      if ( V_f.size() >= 1 )
      {  // compute Newell vector <double>
         Scan_vector normal;
         normal_vector_newell_3(V_f_scan.begin(),V_f_scan.end(),normal);

         // construct and enqueue Nef_polyhedron_3 <Kernel>
         Nef_3 nef (V_f.begin(), V_f.end(), normal, verb);
         if ( !nef.is_empty() )
           {nary_union.add_polyhedron(nef);
            is_nef = true;
         }
      }

      if ( !is_nef )
      {  ++discarded_facets;
         if (verb)
         {  std::cerr << "Hence, discard input facet " << (idx+1)
               << " (enumerated beginning with 1)."
               << " Check semantics!\n" << std::endl;
         }
      }
   }

   nef_union = nary_union.get_union();
   CGAL::Mark_bounded_volumes<Nef_3> mbv (true);
   nef_union.delegate (mbv);

   return discarded_facets;
}


int main(int argc, char *argv[])
{

   if (argc <= 1) {
   Nef_polyhedron N;
   //std::size_t discarded = CGAL::OFF_to_nef_3 (std::cin, N, true);
   std::size_t discarded = HOGE (std::cin, N, true);

   std::cout << "HOGE : " << discarded << std::endl;
   CGAL::draw(N);
   return EXIT_SUCCESS;
   }

  // read OFF file into a polyhedron
  Polyhedron P1, P2;
  std::ifstream ifs1((argc > 1) ? argv[1] : CGAL::data_file_path("meshes/cross_quad.off"));
  ifs1 >> P1;
  std::ifstream ifs2((argc > 2) ? argv[2] : CGAL::data_file_path("meshes/beam.off"));
  ifs2 >> P2;

  // initialize nef from polyhedron
  Nef_polyhedron N1(P1);
  Nef_polyhedron N2(P2);

  Nef_polyhedron NN(N1-N2);
  // draw Nef Polyhedron
  // std::cout << NN << std::endl;

  //CGAL::draw(NN);
  Surface_mesh_ output;
  CGAL::convert_nef_polyhedron_to_polygon_mesh(NN, output);

  //output.vertices();
  std::cout << "ver: " << output.number_of_vertices() << std::endl;
  std::cout << "fcs: " << output.number_of_faces() << std::endl;

  {
      std::cout << "vtx" << std::endl;
      long cntr = 0;
      for(auto it = output.vertices_begin(); it != output.vertices_end(); it++, cntr++) {
          std::size_t ii = (std::size_t)(*it);
          if (ii != cntr) {
              std::cerr << ":failed: ";
          }
          CGAL::Cartesian<double>::Point_3  pt = output.point(*it);
          std::cout << cntr << " " << (*it) << " ";
          std::cout << pt.x() << " ";
          std::cout << pt.y() << " ";
          std::cout << pt.z() << std::endl;
      }
  }
  {
      std::cout << "fc" << std::endl;
      for(auto it = output.faces_begin(); it != output.faces_end(); it++) {
          std::cout << (*it) << " : ";
          auto st = output.halfedge(*it);
          auto nx = output.next(st);
          std::cout << st;
          std::cout << "(" << output.target(st) << ")";
          while(nx != st) {
              std::cout << " " << nx;
              std::cout << "(" << output.target(nx) << ")";
              nx = output.next(nx);
          }
          std::cout << std::endl;
      }
  }

  std::ofstream out;
  out.open("out.off");
  out << output;
  out.close();

  return EXIT_SUCCESS;
}
