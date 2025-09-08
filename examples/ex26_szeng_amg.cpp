#include "mfem.hpp"
#include <iostream>
#include <fstream>

using namespace mfem;
using namespace std;

int main(int argc, char *argv[])
{
    // 1. Initialize MFEM and parse command-line options
    const char *mesh_file = "../data/star.mesh";
    int order = 1;
    bool visualization = true;

    OptionsParser args(argc, argv);
    args.AddOption(&mesh_file, "-m", "--mesh", "Mesh file.");
    args.AddOption(&order, "-o", "--order", "Finite element order.");
    args.AddOption(&visualization, "-vis", "--visualization", "-no-vis", "--no-visualization", "Enable GLVis.");
    args.Parse();
    if (!args.Good())
    {
        args.PrintUsage(cout);
        return 1;
    }
    args.PrintOptions(cout);

    // 2. Load the mesh
    Mesh mesh(mesh_file, 1, 1);
    int dim = mesh.Dimension();

    // 3. Define finite element space
    H1_FECollection fec(order, dim);
    FiniteElementSpace fes(&mesh, &fec);

    cout << "Number of finite element unknowns: " << fes.GetTrueVSize() << endl;

    FunctionCoefficient sigma([](const Vector &x) -> double {
        // Example: Two regions with different conductivities
        return (x(0) < 0.5) ? 10.0 : 1.0;
    });

    // 4. Define bilinear form: A(u,v) = ∫ σ ∇u · ∇v dΩ
    BilinearForm a(&fes);
    a.AddDomainIntegrator(new DiffusionIntegrator(sigma));
    a.Assemble();

    // 5. Define the linear form (RHS)
    LinearForm b(&fes);
    ConstantCoefficient f(1.0);
    b.AddDomainIntegrator(new DomainLFIntegrator(f));
    b.Assemble();

    // 6. Define Dirichlet boundary conditions
    Array<int> ess_bdr(mesh.bdr_attributes.Max());
    ess_bdr = 1; // Apply to all boundaries

    GridFunction x(&fes);
    x = 0.0;

    OperatorPtr A;
    Vector B, X;
    a.FormLinearSystem(ess_bdr, x, b, A, X, B);

        // 7. Define AMG as a preconditioner
    HypreBoomerAMG amg;
    amg.SetPrintLevel(1);  // Print AMG setup
    amg.SetCycleType(HypreBoomerAMG::V_CYCLE);  // Use V-cycle
    amg.SetCoarseningStrategy(HypreBoomerAMG::CLJP);  // Classical coarsening
    amg.SetSmoothingType(1);  // 1 = Jacobi smoothing (can replace with SOR if needed)
    amg.SetMaxLevels(5);  // Limit to 5 multigrid levels
    amg.SetTol(1e-6);

    // 8. Use FGMRES instead of GMRES or PCG
    FGMRESSolver fgmres;
    fgmres.SetPrintLevel(1);
    fgmres.SetMaxIter(200); // Max iterations
    fgmres.SetKDim(30); // Krylov subspace dimension (restart length, mimics GCRO-DR)
    fgmres.SetRelTol(1e-6);
    fgmres.SetPreconditioner(amg);
    fgmres.SetOperator(*A);
    fgmres.Mult(B, X);

        // 9. Recover the solution in finite element space
    a.RecoverFEMSolution(X, b, x);

    // 10. Save and visualize the solution
    ofstream sol_ofs("sol.gf");
    sol_ofs.precision(8);
    x.Save(sol_ofs);

    if (visualization)
    {
        socketstream sol_sock("localhost", 19916);
        sol_sock.precision(8);
        sol_sock << "solution\n" << mesh << x << flush;
    }

    return 0;
}



