#include "mfem.hpp"
#include <fstream>
#include <iostream>
#include <map>
#include <cmath>
#include <tuple>

using namespace std;
using namespace mfem;

// Function to read the conductivity file
std::map<std::tuple<double, double, double>, double> ReadConductivityFile(const std::string &filename)
{
    std::map<std::tuple<double, double, double>, double> conductivity_map;
    std::ifstream infile(filename);
    if (!infile)
    {
        cerr << "Error: Could not open file " << filename << endl;
        exit(1);
    }

    double x, y, z, sigma;
    while (infile >> x >> y >> z >> sigma)
    {
        conductivity_map[std::make_tuple(x, y, z)] = sigma;
    }
    return conductivity_map;
}

// Function for nearest neighbor interpolation of conductivity
double InterpolateSigma(const Vector &x,
                        const std::map<std::tuple<double, double, double>, double> &conductivity_map)
{
    // Get coordinates
    double px = x[0], py = x[1], pz = x[2];

    // Nearest neighbor interpolation
    double min_dist = std::numeric_limits<double>::max();
    double interpolated_sigma = 0.0;
    for (const auto &[coords, sigma] : conductivity_map)
    {
        double cx = std::get<0>(coords);
        double cy = std::get<1>(coords);
        double cz = std::get<2>(coords);
        double dist = sqrt(pow(px - cx, 2) + pow(py - cy, 2) + pow(pz - cz, 2));
        if (dist < min_dist)
        {
            min_dist = dist;
            interpolated_sigma = sigma;
        }
    }

    return interpolated_sigma;
}

// Multigrid solver for the diffusion problem
class DiffusionMultigrid : public GeometricMultigrid
{
private:
    FunctionCoefficient &sigma;

public:
    DiffusionMultigrid(FiniteElementSpaceHierarchy &fespaces, Array<int> &ess_bdr, FunctionCoefficient &sigma_coeff)
        : GeometricMultigrid(fespaces, ess_bdr), sigma(sigma_coeff)
    {
        ConstructCoarseOperatorAndSolver(fespaces.GetFESpaceAtLevel(0));
        for (int level = 1; level < fespaces.GetNumLevels(); ++level)
        {
            ConstructOperatorAndSmoother(fespaces.GetFESpaceAtLevel(level), level);
        }
    }

private:
    void ConstructBilinearForm(FiniteElementSpace &fespace)
    {
        BilinearForm *form = new BilinearForm(&fespace);
        form->SetAssemblyLevel(AssemblyLevel::PARTIAL);
        form->AddDomainIntegrator(new DiffusionIntegrator(sigma));
        form->Assemble();
        bfs.Append(form);
    }

    void ConstructCoarseOperatorAndSolver(FiniteElementSpace &coarse_fespace)
    {
        ConstructBilinearForm(coarse_fespace);

        OperatorPtr opr;
        opr.SetType(Operator::ANY_TYPE);
        bfs[0]->FormSystemMatrix(*essentialTrueDofs[0], opr);
        opr.SetOperatorOwner(false);

        CGSolver *pcg = new CGSolver();
        pcg->SetPrintLevel(-1);
        pcg->SetMaxIter(200);
        pcg->SetRelTol(sqrt(1e-4));
        pcg->SetAbsTol(0.0);
        pcg->SetOperator(*opr.Ptr());

        AddLevel(opr.Ptr(), pcg, true, true);
    }

    void ConstructOperatorAndSmoother(FiniteElementSpace &fespace, int level)
    {
        const Array<int> &ess_tdof_list = *essentialTrueDofs[level];
        ConstructBilinearForm(fespace);

        OperatorPtr opr;
        opr.SetType(Operator::ANY_TYPE);
        bfs[level]->FormSystemMatrix(ess_tdof_list, opr);
        opr.SetOperatorOwner(false);

        Vector diag(fespace.GetTrueVSize());
        bfs[level]->AssembleDiagonal(diag);

        Solver *smoother = new OperatorChebyshevSmoother(*opr, diag, ess_tdof_list, 2);
        AddLevel(opr.Ptr(), smoother, true, true);
    }
};

int main(int argc, char *argv[])
{
    // 1. Parse command-line options
    const char *mesh_file = "../data/star.mesh";
    const char *conductivity_file = "conductivity.txt";
    int geometric_refinements = 0;
    int order_refinements = 2;
    const char *device_config = "cpu";
    bool visualization = true;

    OptionsParser args(argc, argv);
    args.AddOption(&mesh_file, "-m", "--mesh", "Mesh file to use.");
    args.AddOption(&conductivity_file, "-c", "--conductivity-file", "Conductivity file to use.");
    args.AddOption(&geometric_refinements, "-gr", "--geometric-refinements",
                   "Number of geometric refinements.");
    args.AddOption(&order_refinements, "-or", "--order-refinements",
                   "Number of order refinements.");
    args.AddOption(&device_config, "-d", "--device", "Device configuration string.");
    args.AddOption(&visualization, "-vis", "--visualization", "-no-vis",
                   "--no-visualization", "Enable or disable GLVis visualization.");
    args.Parse();
    if (!args.Good())
    {
        args.PrintUsage(cout);
        return 1;
    }
    args.PrintOptions(cout);

    // 2. Enable device
    Device device(device_config);
    device.Print();

    // 3. Read the mesh
    Mesh *mesh = new Mesh(mesh_file, 1, 1);
    int dim = mesh->Dimension();

    // 4. Refine the mesh
    {
        int ref_levels = (int)floor(log(5000. / mesh->GetNE()) / log(2.) / dim);
        for (int l = 0; l < ref_levels; l++)
        {
            mesh->UniformRefinement();
        }
    }

    // 5. Define finite element space hierarchy
    FiniteElementCollection *fec = new H1_FECollection(1, dim);
    FiniteElementSpace *coarse_fespace = new FiniteElementSpace(mesh, fec);
    FiniteElementSpaceHierarchy fespaces(mesh, coarse_fespace, true, true);

    Array<FiniteElementCollection *> collections;
    collections.Append(fec);
    for (int level = 0; level < geometric_refinements; ++level)
    {
        fespaces.AddUniformlyRefinedLevel();
    }
    for (int level = 0; level < order_refinements; ++level)
    {
        collections.Append(new H1_FECollection((int)std::pow(2, level + 1), dim));
        fespaces.AddOrderRefinedLevel(collections.Last());
    }

    // 6. Read conductivity file and create a FunctionCoefficient
    auto conductivity_map = ReadConductivityFile(conductivity_file);
    FunctionCoefficient sigma([&](const Vector &x) -> double {
        return InterpolateSigma(x, conductivity_map);
    });

    // 7. Set up the linear form
    LinearForm *b = new LinearForm(&fespaces.GetFinestFESpace());
    ConstantCoefficient source(1.0); // Example source term
    b->AddDomainIntegrator(new DomainLFIntegrator(source));
    b->Assemble();

    // 8. Define the multigrid solver
    Array<int> ess_bdr(mesh->bdr_attributes.Max());
    ess_bdr = 1; // Set all boundaries to Dirichlet

    DiffusionMultigrid M(fespaces, ess_bdr, sigma);
    M.SetCycleType(Multigrid::CycleType::VCYCLE, 1, 1);

    OperatorPtr A;
    Vector B, X;
    M.FormFineLinearSystem(*b, A, X, B);

    // 9. Solve the system
    PCG(*A, M, B, X, 1, 2000, 1e-12, 0.0);

    // 10. Recover the solution
    M.RecoverFineFEMSolution(X, *b, fespaces.GetFinestFESpace().GetTrueVSize());

    // 11. Save the solution and mesh for visualization
    ofstream mesh_ofs("refined.mesh");
    fespaces.GetFinestFESpace().GetMesh()->Print(mesh_ofs);
    ofstream sol_ofs("sol.gf");
    fespaces.GetFinestFESpace().GetMesh()->Save(sol_ofs);
}
