//                                MFEM
//
// Compile with: make microstrip_issue
//
// This code is about microstrip simulation.
// it compute the electric potential.
// it compute the gradient of the potential which is the electric field.
// It is inspired from MFEM Example 27.
// The mesh is created with gmsh from the file microstrip.geo.
// DL230531.
#include "mfem.hpp"
#include <fstream>
#include <iostream>
#include <miniapps/common/fem_extras.hpp>

using namespace std;
using namespace mfem;

int main(int argc, char *argv[])
{
  // int i;
   const char *mesh_file = "../data/microstrip1_mfem0.mesh";
   int order = 2;
   const char *device_config = "cpu";

   cout << "step #2" << endl;
   // 2. Enable hardware devices such as GPUs, and programming models such as
   //    CUDA, OCCA, RAJA and OpenMP based on command line options.
   Device device(device_config);
   device.Print();

   // 3. Read the mesh from the given mesh file. We can handle triangular,
   //    quadrilateral, tetrahedral, hexahedral, surface and volume meshes with
   //    the same code.
   Mesh mesh(mesh_file, 1, 1);
   mesh.UniformRefinement();
   int dim = mesh.Dimension();

  // 3. Define a finite element space on the mesh. Here we use 
    //    continuous Lagrange finite elements.
    FiniteElementCollection *fec = (FiniteElementCollection*)new H1_FECollection(order, dim);
    FiniteElementSpace fespace(&mesh, fec);
    int size = fespace.GetTrueVSize();
    cout << "Number of finite element unknowns: " << size << endl;
   
   // 4. Create "marker arrays" to define the portions of boundary associated
   //    with each type of boundary condition. These arrays have an entry
   //    corresponding to each boundary attribute.  Placing a '1' in entry i
   //    marks attribute i+1 as being active, '0' is inactive.
   //    in this case there are omly dirichelet boundary.
   
   Array<int> dbc_bdr(mesh.bdr_attributes.Max());
   assert(mesh.bdr_attributes.Max()==2);
   dbc_bdr = 0; dbc_bdr[0] = 1; dbc_bdr[1] = 1;

   Array<int> ess_tdof_list(0);
   if (mesh.bdr_attributes.Size())
   {
      // For a continuous basis the linear system must be modified to enforce an
      // essential (Dirichlet) boundary condition. 
      fespace.GetEssentialTrueDofs(dbc_bdr, ess_tdof_list);
   }
  // pcb dielectric 5, copper 1 and air 1.
   double CoeffArray[]={0.0, 0.0, 5.0, 1.0, 1.0};
   Vector CoeffVector(CoeffArray, 5);
   PWConstCoefficient Coeff(CoeffVector);

   cout << "step #6" << endl;
   // 6. Define the solution vector u as a finite element grid function
   //    corresponding to fespace. Initialize u with initial guess of zero.
   GridFunction u(&fespace);
   u = 0.0;

   // 7. Set up the bilinear form a(.,.) on the finite element space
   //    corresponding to the Laplacian operator -Delta, by adding the Diffusion
   //    domain integrator.
   BilinearForm a(&fespace);
   a.AddDomainIntegrator(new DiffusionIntegrator(Coeff));
   a.Assemble();

   // 8. Assemble the linear form for the right hand side vector.
   LinearForm b(&fespace);

   // Set the Dirichlet values in the solution vector
   // the trace is at 1V, the contour is 0V.
   double BoundaryCoeffArray[]={0.0, 1.0};
   Vector BoundaryCoeffVector(BoundaryCoeffArray, 2);
   PWConstCoefficient BoundaryCoeff(BoundaryCoeffVector);
   u.ProjectBdrCoefficient(BoundaryCoeff, dbc_bdr);
   b.Assemble();

   cout << "step #9" << endl;
   // 9. Construct the linear system.
   {
     OperatorPtr A;
     Vector B, X;
     a.FormLinearSystem(ess_tdof_list, u, b, A, X, B);
   
     // 10. Define a simple symmetric Gauss-Seidel preconditioner and use it to
     //     solve the system AX=B with PCG in the symmetric case, and GMRES in the
     //     non-symmetric one.
     GSSmoother M((SparseMatrix&)(*A));
     PCG(*A, M, B, X, 1, 500, 1e-12, 0.0);
     cout << "step #12" << endl;
     // 12. Recover the grid function corresponding to U. This is the local finite
     //     element solution.
     a.RecoverFEMSolution(X, b, u);
   }

   // 14. Save the refined mesh and the solution. This output can be viewed
   //     later using GLVis: "glvis -m refined.mesh -g sol.gf".
   {
      ofstream mesh_ofs("refined.mesh");
      mesh_ofs.precision(8);
      mesh.Print(mesh_ofs);
      ofstream sol_ofs("sol.gf");
      sol_ofs.precision(8);
      u.Save(sol_ofs);
   }

   // 15. Send the potential solution by socket to a GLVis server.
   string title_str = "H1";
   char vishost[] = "localhost";
   int  visport   = 19916;
   socketstream sol_sock(vishost, visport);
   sol_sock.precision(8);
   sol_sock << "solution\n" << mesh << u
	    << "window_title '" << title_str << " Solution'"
	    << " keys 'mmc'" << flush;
   
   cout << "step compute gradient" << endl;
   // This section computes the gradient (field).
   ND_FECollection nd_fec(order, dim);
   FiniteElementSpace fespace_nd(&mesh, &nd_fec);
   GridFunction dT(&fespace_nd);
   {
     DiscreteLinearOperator grad(&fespace, &fespace_nd);
     grad.AddDomainInterpolator(new GradientInterpolator);
     grad.Assemble();
     grad.Finalize();
     grad.Mult(u, dT);
   }
   
   RT_FECollection rt_fec(order, dim);
   FiniteElementSpace fespace_rt(&mesh, &rt_fec);
   GridFunction D(&fespace_rt);
   {
     LinearForm epsdT(&fespace_rt);
     MixedBilinearForm epsGrad(&fespace, &fespace_rt);
     epsGrad.AddDomainIntegrator(new MixedVectorGradientIntegrator(Coeff));
     epsGrad.Assemble();
     epsGrad.Finalize();
     epsGrad.Mult(u, epsdT);

     BilinearForm m_rt(&fespace_rt);
     m_rt.AddDomainIntegrator(new VectorFEMassIntegrator);
     m_rt.Assemble();
     m_rt.Finalize();

     Array<int> ess_tdof_rt_list;
     OperatorPtr A;
     Vector B, X;

     D = 0.0;
     m_rt.FormLinearSystem(ess_tdof_rt_list, D, epsdT, A, X, B);

     GSSmoother M((SparseMatrix&)(*A));
     PCG(*A, M, B, X, 1, 500, 1e-12, 0.0);
     m_rt.RecoverFEMSolution(X, epsdT, D);
     D *= -1.0;
   }
   
   //Send the gradient solution by socket to a GLVis server.
   {
     string title_str = "Gradient, Field";
     char vishost[] = "localhost";
     int  visport   = 19916;
     socketstream sol_sock(vishost, visport);
     sol_sock.precision(8);
     sol_sock << "solution\n" << mesh << dT
	      << "window_title '" << title_str << " Solution'"
	      << " keys 'mmcvv'" << flush;
   }
   {
     string title_str = "Displacement Field";
     char vishost[] = "localhost";
     int  visport   = 19916;
     socketstream sol_sock(vishost, visport);
     sol_sock.precision(8);
     sol_sock << "solution\n" << mesh << D
	      << "window_title '" << title_str << " Solution'"
	      << " keys 'mmcvv'" << flush;
   }

   // from ex5.cpp  
   // 14. Save data in the VisIt format
   VisItDataCollection visit_dc("microstrip", &mesh);
   visit_dc.RegisterField("potential", &u);
   visit_dc.RegisterField("gradient", &dT);
   visit_dc.Save();

   cout << "step compute integral" << endl;
   //Compute the integral of the gradient on the boundary,
   //this will be the charge.
   // With the help of https://github.com/mfem/mfem/issues/993
   // and from from volta_solver.cpp
   {
     LinearForm rt_surf_int_(&fespace_rt);
     Array<int> bdr_marker(2); bdr_marker[0]=1; bdr_marker[1]=0;
     rt_surf_int_.AddBoundaryIntegrator(new VectorFEBoundaryFluxLFIntegrator,
					bdr_marker);
     rt_surf_int_.Assemble();
     double charge_D1 = rt_surf_int_(D);
     cout << endl<<"charge_D1 = "<<charge_D1<<endl;
   }
   {
     LinearForm rt_surf_int_(&fespace_rt);
     Array<int> bdr_marker(2); bdr_marker[0]=0; bdr_marker[1]=1;
     rt_surf_int_.AddBoundaryIntegrator(new VectorFEBoundaryFluxLFIntegrator,
					bdr_marker);
     rt_surf_int_.Assemble();
     double charge_D2 = rt_surf_int_(D);
     cout << endl<<"charge_D2 = "<<charge_D2<<endl;
   }

   // 16. Free the used memory.
   delete fec;
   mesh.Save("microstrip1_mfem1.mesh");
   return 0;
}
