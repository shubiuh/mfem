
#include <mfem.hpp>
#include <miniapps/common/fem_extras.hpp>
#include <fstream>
#include <iostream>
#include <string.h>

/* These definitions fit with the mesh file ProxRoundWires3d.msh */
#define air 1 //domain
#define conductorright 2 //domain
#define conductorleft 3  //domain
#define airsurffront 4  //then all boundary
#define airsurfback 5
#define conductorleftsurffront 6
#define conductorleftsurfback 7
#define conductorrightsurffront 8
#define conductorrightsurfback 9
#define airsurfaround 10

using namespace std;
using namespace mfem;

class ScalarPotential
{
   protected:
      Mesh *mesh;      // pointer to mfem mesh file.
      const char *mesh_file;   // mes file name.
      int order;       // order for potential elements.
      int refineTo;    // Default 1 cause no refinement.
      int RefineCount; // Number of time refine is called.
      int dim;         // mesh dimension.

      FiniteElementCollection *fec;
      FiniteElementSpace *fes;
      int size; // number of vector true (conforming) dofs.

      Array<int> *dbc_bdr;   //dirichelet boundary marker array.
      Array<int> *ess_tdof_list;  //essential true degree of freedom list.
      int ess_tdof_list_size;   
      PWConstCoefficient *BoundaryCoeff;
      PWConstCoefficient *DomainCoeff;
      GridFunction *gf;
      BilinearForm *blf;
      LinearForm *lf;
      OperatorPtr A;  
      Vector *B, *X;
public:
      ScalarPotential();
      void ReadArgs(int argc, char *(argv[]));
      void LoadMesh();
      void CreateFESpace();
      void PrepareDomainCoeff();
      void PrepareDiricheletBoundaryCond();
      void CreateFunctional();
      void Solve();
      void SaveAndDisplay();
      ~ScalarPotential();
    
    Mesh *GetMeshPtr() {return mesh;};
    FiniteElementSpace *GetFes() {return fes;};
    PWConstCoefficient *GetDomainCoeff() {return DomainCoeff;};
    GridFunction *GetGridFunction() {return gf;};
    int GetOrder() {return order;};
    
};


class EFieldFromVoltagePotential
{
   protected:
      Mesh *mesh;      // pointer to mfem mesh file.
      FiniteElementSpace *pfes;
      GridFunction *pgf;

      int rt_order;       // order for RT vector elements.
      int dim;         // mesh dimension.

      FiniteElementCollection *fec;
      FiniteElementSpace *fes;
      int size; // number of vector true (conforming) dofs.
      Array<int> *ess_tdof_list;  //essential true degree of freedom list.
      int ess_tdof_list_size;   
      PWConstCoefficient *DomainCoeff;
      GridFunction *gf;
      MixedBilinearForm *mblf;
      BilinearForm *blf;
      LinearForm *lf;
      OperatorPtr A;  
      Vector *B, *X;
public:
      EFieldFromVoltagePotential(GridFunction *pgf, PWConstCoefficient *PWCC);
      void LoadMesh();
      void CreateFESpace();
      void PrepareDiricheletBoundaryCond();
      void CreateFunctional();
      void Solve();
      void SaveAndDisplay();
      ~EFieldFromVoltagePotential();
    
};

  ScalarPotential::ScalarPotential()
  {
     mesh_file = "../data/proxroundwires3d.msh"; //default mesh file.
     order = 2; // default order for potential elements.
     refineTo = 1; // Default 1 cause no refinement.

     //  Enable hardware devices such as GPUs, and programming models such as
     //    CUDA, OCCA, RAJA and OpenMP based on command line options.
     Device device("cpu");
     device.Print();
  }

  void ScalarPotential::ReadArgs(int argc, char *(argv[]))
  {
     OptionsParser args(argc, argv);
      args.AddOption(&mesh_file, "-m", "--mesh", "Mesh file to use.");
      args.AddOption(&order, "-o", "--order", "Finite element polynomial degree");
      args.AddOption(&refineTo, "-rt", "--refineto", "Refine to _ elements");
      args.ParseCheck();  
  }

void ScalarPotential::LoadMesh()
{
   //  Read the mesh from the given mesh file. We can handle triangular,
   //  quadrilateral, tetrahedral, hexahedral, surface and volume meshes with
   //  the same code.
   mesh = new Mesh(mesh_file, 1, 1);

   //Refine the mesh, better to refine with the generator tool to improve section borders.
   RefineCount = 0;
   while(mesh->GetNE()<refineTo) {
      mesh->UniformRefinement();
      RefineCount++;
   } 
   cout << "Refine " << RefineCount << " times."<<endl; 
   
   dim = mesh->Dimension();

   cout << "mesh.Dimension() = "<< mesh->Dimension() << endl;
   cout << "mesh.GetNE() = "<< mesh->GetNE() << endl;
   cout << "mesh.GetNBE() = "<< mesh->GetNBE() << endl;
   cout << "mesh.GetNEdges() = "<< mesh->GetNEdges() << endl;
   cout << "mesh.GetNFaces() = "<< mesh->GetNFaces() << endl;
   cout << "mesh.bdr_attributes.Max() = "<< mesh->bdr_attributes.Max() << endl;

}

void ScalarPotential::CreateFESpace()
{
   //  Define a finite element space on the mesh. Here we use 
   //  continuous Lagrange finite elements, since we seek a scalar.
    fec = (FiniteElementCollection*)new H1_FECollection(order, dim);
    fes = new FiniteElementSpace(mesh, fec);
    int size = fes->GetTrueVSize();
    cout << "Number of finite element unknowns: " << size << endl;

   //  Define the solution vector u as a finite element grid function
   //  corresponding to fes. Initialize u with initial guess of zero.
   gf = new GridFunction(fes);
   *gf = 0.0;
}

void ScalarPotential::PrepareDiricheletBoundaryCond()
{
   // 4. Create "marker arrays" to define the portions of boundary associated
   //    with each type of boundary condition. These arrays have an entry
   //    corresponding to each boundary attribute.  Placing a '1' in entry i
   //    marks attribute i+1 as being active, '0' is inactive.
   //    in this case there are only dirichelet boundary.
      
   dbc_bdr = new Array<int>(mesh->bdr_attributes.Max());
   cout << "mesh.bdr_attributes.Max() = " << mesh->bdr_attributes.Max() << endl;
   assert(mesh->bdr_attributes.Max()==6);  // to fit the mesh file.
   *dbc_bdr = 1;
   (*dbc_bdr)[conductorright-1]=0;
   (*dbc_bdr)[conductorleft-1]=0;
   (*dbc_bdr)[air-1]=0;
   // both airsurffront and airsurfback should not be forced boundary value
   // to avoid very large gradient at the boundary.
   (*dbc_bdr)[airsurffront-1]=0; 
   (*dbc_bdr)[airsurfback-1]=0;
   if(0) dbc_bdr->Print(cout);

   ess_tdof_list = new Array<int>(0);
   if (mesh->bdr_attributes.Size())
   {
      // For a continuous basis the linear system must be modified to enforce an
      // essential (Dirichlet) boundary condition. 
      fes->GetEssentialTrueDofs(*dbc_bdr, *ess_tdof_list);
   }

   ess_tdof_list_size = ess_tdof_list->Size();
   cout << "ess_tdof_list_size = " << ess_tdof_list_size << endl;

   // Set the Dirichlet values in the solution vector
   double *CoeffArray;
   CoeffArray = new double[mesh->bdr_attributes.Max()];
   CoeffArray[air-1]=0.0;
   CoeffArray[conductorright-1]=0.0;
   CoeffArray[conductorleft-1]=0.0;
   CoeffArray[airsurffront-1]=0.0;
   CoeffArray[airsurfback-1]=0.0;
   CoeffArray[conductorleftsurffront-1]=-1.0;  // -1.0 Volt.
   CoeffArray[conductorleftsurfback-1]=1.0;    // +1.0 Volt.
   CoeffArray[conductorrightsurffront-1]=1.0;
   CoeffArray[conductorrightsurfback-1]=-1.0;
   CoeffArray[airsurfaround-1]=0.0;

   Vector BoundaryCoeffVector(CoeffArray, mesh->bdr_attributes.Max());
   delete CoeffArray;
   BoundaryCoeff = new PWConstCoefficient(BoundaryCoeffVector);
   gf->ProjectBdrCoefficient(*BoundaryCoeff, *dbc_bdr);
}

void ScalarPotential::PrepareDomainCoeff()
{
  //conduction air, copper, copper.
   double CoeffArray[] = {0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
   assert((sizeof(CoeffArray)/sizeof(CoeffArray[0])) == mesh->bdr_attributes.Max());
   Vector CoeffVector(CoeffArray, mesh->bdr_attributes.Max());
   DomainCoeff = new PWConstCoefficient(CoeffVector);

}

void ScalarPotential::CreateFunctional()
{
   //  Set up the bilinear form a(.,.) on the finite element space
   //    corresponding to the Laplacian operator -Delta, by adding the Diffusion
   //    domain integrator.
   blf = new BilinearForm(fes);
   blf->AddDomainIntegrator(new DiffusionIntegrator(*DomainCoeff));
   blf->Assemble();

   //  Assemble the linear form for the right hand side vector.
   lf = new LinearForm(fes);
   lf->Assemble();

   //  Construct the linear system.
   //A = new OperatorPtr;
   B = new Vector;
   X = new Vector;
   blf->FormLinearSystem(*ess_tdof_list, *gf, *lf, A, *X, *B);

}

void ScalarPotential::Solve()
{
      //  Construct the linear system.
   //A = new OperatorPtr;
   B = new Vector;
   X = new Vector;
   blf->FormLinearSystem(*ess_tdof_list, *gf, *lf, A, *X, *B);
   
   //  Define a simple symmetric Gauss-Seidel preconditioner and use it to
   //     solve the system AX=B with PCG in the symmetric case, and GMRES in the
   //     non-symmetric one.
        GSSmoother M((SparseMatrix&)(*A));
         PCG(*A, M, *B, *X, 1, 500, 1e-12, 0.0);
 
//  Recover the grid function corresponding to U. This is the local finite
//     element solution.
   blf->RecoverFEMSolution(*X, *lf, *gf);
}

void ScalarPotential::SaveAndDisplay()
{
   
   // 15. Send the potential solution by socket to a GLVis server.
      string title_str = "H1";
      char vishost[] = "localhost";
      int  visport   = 19916;
      socketstream sol_sock(vishost, visport);
      sol_sock.precision(8);
      sol_sock << "solution\n" << *mesh << *gf
               << "window_title '" << title_str << " Solution'"
               << " keys 'mmc'" << flush;

// 14. Save data in the VisIt format
   VisItDataCollection visit_dc("two_wires_poy", mesh);
   visit_dc.RegisterField("potential", gf);
   visit_dc.Save();

   // 15. Save data in the ParaView format
   ParaViewDataCollection paraview_dc("two_wires_pw_pot", mesh);
   paraview_dc.SetPrefixPath("ParaView");
   paraview_dc.SetLevelsOfDetail(order);
   paraview_dc.SetCycle(0);
   paraview_dc.SetDataFormat(VTKFormat::BINARY);
   paraview_dc.SetHighOrderOutput(true);
   paraview_dc.SetTime(0.0); // set the time
   paraview_dc.RegisterField("potential", gf);
   paraview_dc.Save();


}

ScalarPotential::~ScalarPotential()
{
    
    
   delete dbc_bdr;
   
   delete ess_tdof_list;

   delete DomainCoeff;

   delete gf;
   delete blf;

   delete lf;

   delete BoundaryCoeff;
  
   delete B;
   delete X;


   delete fec;
   delete fes;
   delete mesh;
}


  EFieldFromVoltagePotential::EFieldFromVoltagePotential(GridFunction *pgf, PWConstCoefficient *PWCC)
  {
     pfes = pgf->FESpace();
     rt_order = pfes->FEColl()->GetOrder()-1; 
     mesh = pfes->GetMesh();
     DomainCoeff = PWCC;
     //  Enable hardware devices such as GPUs, and programming models such as
     //    CUDA, OCCA, RAJA and OpenMP based on command line options.
     Device device("cpu");
     device.Print();
  }

void EFieldFromVoltagePotential::LoadMesh()
{
   //  the mesh pointer is already set by the constructor.
   
   dim = mesh->Dimension();


}

void EFieldFromVoltagePotential::CreateFESpace()
{
// This section computes the gradient (field) using raviart thomas basis function.
   fec = new RT_FECollection(rt_order, dim);
   fes = new FiniteElementSpace(mesh, fec);
   gf = new GridFunction(fes);  //trial space is RT.
   
     lf = new LinearForm(fes);
     mblf = new MixedBilinearForm(pfes, fes);
     mblf->AddDomainIntegrator(new MixedVectorGradientIntegrator(/* *DomainCoeff */));
     mblf->Assemble();
     mblf->Finalize();
     cout << pgf->Size() << endl;
     mblf->Mult(*pgf, *lf);

     blf = new BilinearForm(fes);
     blf->AddDomainIntegrator(new VectorFEMassIntegrator);
     blf->Assemble();
     blf->Finalize();


     Array<int> ess_tdof_rt_list;
     B = new Vector;
     X = new Vector;
     

     *gf = 0.0;
     blf->FormLinearSystem(ess_tdof_rt_list, *gf, *lf, A, *X, *B);


    int size = fes->GetTrueVSize();
    cout << "Number of finite element unknowns: " << size << endl;

 }

void EFieldFromVoltagePotential::Solve()
{
     GSSmoother M((SparseMatrix&)(*A));
     PCG(*A, M, *B, *X, 1, 500, 1e-12, 0.0);
     blf->RecoverFEMSolution(*X, *lf, *gf);
     *gf *= -1.0;
}

void EFieldFromVoltagePotential::SaveAndDisplay()
{
    
     string title_str = "Electric Field";
     char vishost[] = "localhost";
     int  visport   = 19916;
     socketstream sol_sock(vishost, visport);
     sol_sock.precision(8);
     sol_sock << "solution\n" << *mesh << *gf
	      << "window_title '" << title_str << " Solution'"
	      << " keys 'mmcvv'" << flush;
   

// from ex5.cpp  
// 14. Save data in the VisIt format
   VisItDataCollection visit_dc("two_wires_Efield", mesh);
   visit_dc.RegisterField("Efield", gf);
   visit_dc.Save();

   // 15. Save data in the ParaView format
   ParaViewDataCollection paraview_dc("two_wires_pw_Efield", mesh);
   paraview_dc.SetPrefixPath("ParaView");
   paraview_dc.SetLevelsOfDetail(rt_order);
   paraview_dc.SetCycle(0);
   paraview_dc.SetDataFormat(VTKFormat::BINARY);
   paraview_dc.SetHighOrderOutput(true);
   paraview_dc.SetTime(0.0); // set the time
   paraview_dc.RegisterField("Efield",gf);
   paraview_dc.Save();

}


int main(int argc, char *argv[])
{
   ScalarPotential *SP = new ScalarPotential();
   SP->ReadArgs(argc, argv);
   SP->LoadMesh();
   SP->CreateFESpace();
   SP->PrepareDiricheletBoundaryCond();
   SP->PrepareDomainCoeff();
   SP->CreateFunctional();
   SP->Solve();
   SP->SaveAndDisplay();

   EFieldFromVoltagePotential *EF = new   EFieldFromVoltagePotential(SP->GetGridFunction(), SP->GetDomainCoeff());
   EF->LoadMesh();
   EF->CreateFESpace();  
   EF->Solve();
   EF->SaveAndDisplay();



   delete SP;
}