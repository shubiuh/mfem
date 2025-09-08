//                                MFEM Example 0
//
// Compile with: make ex0
//
// Sample runs:  ex0
//               ex0 -m ../data/fichera.mesh
//               ex0 -m ../data/square-disc.mesh -o 2
//
// Description: This example code demonstrates the most basic usage of MFEM to
//              define a simple finite element discretization of the Laplace
//              problem -Delta u = 1 with zero Dirichlet boundary conditions.
//              General 2D/3D mesh files and finite element polynomial degrees
//              can be specified by command line options.

#include "mfem.hpp"
#include "nanoflann.hpp"
#include <vector>
#include <array>
#include <fstream>
#include <iostream>

using namespace std;
using namespace mfem;

// 3D grid-based storage for conductivity values
struct ConductivityGrid {
    int nx, ny, nz;  // Grid dimensions
    double x_min, x_max, y_min, y_max, z_min, z_max;  // Domain bounds
    vector<double> values;  // Flattened 3D array for conductivity values

    // Add a default constructor
    ConductivityGrid() : nx(0), ny(0), nz(0),
                        x_min(0), x_max(0),
                        y_min(0), y_max(0),
                        z_min(0), z_max(0) {}

    ConductivityGrid(int nx_, int ny_, int nz_,
                     double x_min_, double x_max_,
                     double y_min_, double y_max_,
                     double z_min_, double z_max_)
        : nx(nx_), ny(ny_), nz(nz_),
          x_min(x_min_), x_max(x_max_),
          y_min(y_min_), y_max(y_max_),
          z_min(z_min_), z_max(z_max_)
    {
        values.resize(nx * ny * nz, 1.0);  // Default sigma = 1.0
    }

    void update(int nx_, int ny_, int nz_,
               double x_min_, double x_max_,
               double y_min_, double y_max_,
               double z_min_, double z_max_,
               std::vector<double> sigma_values_in) {
        nx = nx_;
        ny = ny_;
        nz = nz_;
        x_min = x_min_;
        x_max = x_max_;
        y_min = y_min_;
        y_max = y_max_;
        z_min = z_min_;
        z_max = z_max_;
        values.resize(nx * ny * nz, 1.0);  // Default sigma = 1.0
        values = sigma_values_in;
    }

    void TransformValuse(double multiplier, double offset) {
        for (auto& value : values) {
            value = value * multiplier + offset;
            //if (value > 30) {
            //    value = 10;
            //}
            //else
            //{
            //    value = 0;
            //}
        }
    }

    // Convert 3D indices to 1D array index
    inline int index(int i, int j, int k) const {
        // return i + nx * (j + ny * k);
        return k + nz * (j + ny * i);
    }

    // Load conductivity data from a file
    void LoadFromFile(const string &filename) {
        ifstream infile(filename);
        if (!infile) {
            cerr << "Error: Could not open " << filename << endl;
            exit(1);
        }

        double x, y, z, sigma;
        while (infile >> x >> y >> z >> sigma) {
            int i = static_cast<int>((x - x_min) / (x_max - x_min) * (nx - 1));
            int j = static_cast<int>((y - y_min) / (y_max - y_min) * (ny - 1));
            int k = static_cast<int>((z - z_min) / (z_max - z_min) * (nz - 1));

            i = std::min(std::max(i, 0), nx - 1);
            j = std::min(std::max(j, 0), ny - 1);
            k = std::min(std::max(k, 0), nz - 1);

            values[index(i, j, k)] = sigma*100+10;
        }
    }

    // Trilinear interpolation of sigma(x, y, z)
    double Interpolate(double x, double y, double z, std::string interplationscheme = "nearest") const {
        // Normalize x, y, z to grid indices
        double fx = (x - x_min) / (x_max - x_min) * (nx - 1);
        double fy = (y - y_min) / (y_max - y_min) * (ny - 1);
        double fz = (z - z_min) / (z_max - z_min) * (nz - 1);

        // Find bounding indices
        int i0 = static_cast<int>(std::round(fx)), i1 = i0 + 1;
        int j0 = static_cast<int>(std::round(fy)), j1 = j0 + 1;
        int k0 = static_cast<int>(std::round(fz)), k1 = k0 + 1;

        // Clamp indices to valid range
        i0 = std::max(i0, 0);
        j0 = std::max(j0, 0);
        k0 = std::max(k0, 0);
        i1 = std::min(i1, nx - 1);
        j1 = std::min(j1, ny - 1);
        k1 = std::min(k1, nz - 1);

        // Get interpolation weights
        double wx = fx - i0, wy = fy - j0, wz = fz - k0;

        // Trilinear interpolation
        if (interplationscheme == "nearest")
        {
            return values[index(i0, j0, k0)];
        } else if (interplationscheme == "linear")
        {
            double c000 = values[index(i0, j0, k0)];
            double c100 = values[index(i1, j0, k0)];
            double c010 = values[index(i0, j1, k0)];
            double c110 = values[index(i1, j1, k0)];
            double c001 = values[index(i0, j0, k1)];
            double c101 = values[index(i1, j0, k1)];
            double c011 = values[index(i0, j1, k1)];
            double c111 = values[index(i1, j1, k1)];

            return (1 - wx) * (1 - wy) * (1 - wz) * c000 +
                wx * (1 - wy) * (1 - wz) * c100 +
                (1 - wx) * wy * (1 - wz) * c010 +
                wx * wy * (1 - wz) * c110 +
                (1 - wx) * (1 - wy) * wz * c001 +
                wx * (1 - wy) * wz * c101 +
                (1 - wx) * wy * wz * c011 +
                wx * wy * wz * c111;
        }
    }
};

void ReadCOMSOLGrid(const string& filename, ConductivityGrid& grid)
{
    ifstream infile(filename);
    if (!infile.is_open())
    {
        throw runtime_error("Error: Could not open file " + filename);
    }

    string line;
    std::vector<double> x_coords, y_coords, z_coords;
    std::vector<double> sigma_values;

    enum Section
    {
        NONE,
        GRID,
        DATA
    } section = NONE;

    while (getline(infile, line))
    {
        // Determine section
        if (line.find("% Grid") != string::npos)
        {
            section = GRID;
            continue;
        }
        else if (line.find("% Data") != string::npos)
        {
            section = DATA;
            continue;
        }

        // Parse grid or data values
        if (section == GRID)
        {
            istringstream iss(line);
            std::vector<double> coords((istream_iterator<double>(iss)), istream_iterator<double>());
            if (x_coords.empty())
                x_coords = coords;
            else if (y_coords.empty())
                y_coords = coords;
            else
                z_coords = coords;
        }
        else if (section == DATA)
        {
            istringstream iss(line);
            std::vector<double> row((istream_iterator<double>(iss)), istream_iterator<double>());
            sigma_values.insert(sigma_values.end(), row.begin(), row.end());
        }
    }

    // Verify that the data matches the grid dimensions
    int nx = x_coords.size();
    int ny = y_coords.size();
    int nz = z_coords.size();

    if (sigma_values.size() != nx * ny * nz)
    {
        throw runtime_error("Error: Mismatch between grid dimensions and data size.");
    }

    // Initialize the ConductivityGrid structure
    grid.update(nx, ny, nz,
        x_coords.front(), x_coords.back(),
        y_coords.front(), y_coords.back(),
        z_coords.front(), z_coords.back(),
        sigma_values);
    grid.TransformValuse(1, 0);
}


int main(int argc, char *argv[])
{
   StopWatch sw;
   sw.Start();
   // 1. Parse command line options.
   string mesh_file = "../data/inline-hex-51.mesh";
   int order = 1;
   bool visualization = true;
   bool paraview_output = false;

   OptionsParser args(argc, argv);
   args.AddOption(&mesh_file, "-m", "--mesh", "Mesh file to use.");
   args.AddOption(&order, "-o", "--order", "Finite element polynomial degree");
   args.AddOption(&visualization, "-vis", "--visualization", "-no-vis",
               "--no-visualization",
               "Enable or disable GLVis visualization.");
   args.AddOption(&paraview_output, "-pv", "--paraview", "-no-pv",
               "--no-paraview",
               "Enable or disable Paraview output.");
   args.ParseCheck();

   // 2. Read the mesh from the given mesh file, and refine once uniformly.
   Mesh mesh(mesh_file);
   //mesh.UniformRefinement();
   std::cout << "mesh->attributes.Max(): " << mesh.attributes.Max() << std::endl;
   std::cout << "mesh->bdr_attributes.Max(): " << mesh.bdr_attributes.Max() << std::endl;

   // 3. Define a finite element space on the mesh. Here we use H1 continuous
   //    high-order Lagrange finite elements of the given order.
   int dim = mesh.Dimension();
   H1_FECollection fec(order, dim);
   FiniteElementSpace fespace(&mesh, &fec);
   cout << "Number of unknowns: " << fespace.GetTrueVSize() << endl;

   // 4. Extract the list of all the boundary DOFs. These will be marked as
   //    Dirichlet in order to enforce zero boundary conditions.
   Array<int> ess_tdof_list, ess_bdr(mesh.bdr_attributes.Max());
   ess_bdr = 0;
   int top_marker = 5;
   int bottom_marker = 0;

   // 5. Define the solution x as a finite element grid function in fespace. Set
   //    the initial guess to zero, which also sets the boundary conditions.
   GridFunction x(&fespace);
   x = 0.0;
   ess_bdr[top_marker] = 1; // set top to be essential
   ess_bdr[bottom_marker] = 0; // set bottom to be essential
   ConstantCoefficient one(1.0);
   x.ProjectBdrCoefficient(one, ess_bdr);
   ess_bdr[top_marker] = 0; // set top to be essential
   ess_bdr[bottom_marker] = 1; // set bottom to be essential
   ConstantCoefficient zero(0.0);
   x.ProjectBdrCoefficient(zero, ess_bdr);
   ess_bdr[top_marker] = 1; // set top to be essential
   ess_bdr[bottom_marker] = 1; // set bottom to be essential
   fespace.GetEssentialTrueDofs(ess_bdr, ess_tdof_list);

   // Create a FunctionCoefficient using the interpolation function
   ConductivityGrid conductivity_grid = ConductivityGrid();
   //ReadCOMSOLGrid("data51.txt", conductivity_grid);
   ReadCOMSOLGrid("case51_new.txt", conductivity_grid);

   //conductivity_grid.LoadFromFile("conductivity.txt");

   FunctionCoefficient sigma([&](const Vector& x) -> double {
       return conductivity_grid.Interpolate(x[0], x[1], x[2]);
       });

   // 6. Set up the linear form b(.) corresponding to the right-hand side.
//    ConstantCoefficient one(1.0);
//    ConstantCoefficient zero(0.0);
   LinearForm b(&fespace);
   b.AddDomainIntegrator(new DomainLFIntegrator(zero));
//    b.AddDomainIntegrator(new DomainLFIntegrator(sigma));
   b.Assemble();

   // Create and save a grid function for sigma
   GridFunction sigma_gf(&fespace);
   sigma_gf.ProjectCoefficient(sigma);
   sigma_gf.Save("sigma.gf");
   

   // 7. Set up the bilinear form a(.,.) corresponding to the -Delta operator.
   BilinearForm a(&fespace);
   a.AddDomainIntegrator(new DiffusionIntegrator(sigma));
   a.Assemble();

   // 8. Form the linear system A X = B. This includes eliminating boundary
   //    conditions, applying AMR constraints, and other transformations.
   SparseMatrix A;
   Vector B, X;
   a.FormLinearSystem(ess_tdof_list, x, b, A, X, B);

   // 9. Solve the system using PCG with symmetric Gauss-Seidel preconditioner.
   GSSmoother M(A);
   PCG(A, M, B, X, 1, 500, 1e-16, 0.0);

   // 10. Recover the solution x as a grid function and save to file. The output
   //     can be viewed using GLVis as follows: "glvis -m mesh.mesh -g sol.gf"
   a.RecoverFEMSolution(X, b, x);
   x.Save("sol.gf");
   mesh.Save("mesh.mesh");

   // Create and visualize boundary attributes as a grid function
   /*GridFunction bdr_attr(&fespace);
   bdr_attr = 0.0;
   for (int i = 0; i < mesh.GetNBE(); i++)
   {
      Array<int> dofs;
      mesh.GetBdrElementVertices(i, dofs);
      for (int j = 0; j < dofs.Size(); j++)
      {
         if (dofs[j] >= 0)
         {
            bdr_attr(dofs[j]) = mesh.GetBdrAttribute(i);
         }
      }
   }*/
   // Save boundary attributes visualization
   //bdr_attr.Save("bdr_attr.gf");

   // Visualize boundary attributes
   /*if (visualization)
   {
      socketstream bdr_sock("localhost", 19916);
      bdr_sock.precision(8);
      bdr_sock << "solution\n" << mesh << bdr_attr << flush;
   }*/

   sw.Stop();
   double elapsed_time = sw.RealTime();
   std::cout << "Elapsed time: " << elapsed_time << " seconds." << std::endl;

   // 14. Send the solution by socket to a GLVis server.
    if (visualization)
    {
        char vishost[] = "localhost";
        int  visport   = 19916;
        socketstream sol_sock(vishost, visport);
        sol_sock.precision(8);
        sol_sock << "solution\n" << mesh << x << flush;
    }

    ParaViewDataCollection* pd = NULL;
    if (paraview_output)
    {
        pd = new ParaViewDataCollection("Example0", &mesh);
        pd->SetPrefixPath("ParaView");
        pd->RegisterField("solution", &x);
        pd->SetLevelsOfDetail(order);
        pd->SetDataFormat(VTKFormat::BINARY);
        pd->SetHighOrderOutput(true);
        pd->SetCycle(0);
        pd->SetTime(0.0);
        pd->Save();
    }

    // 17. compute gradient (field)
    std::cout << "step compute gradient" << std::endl;
    ND_FECollection nd_fec(order, dim);
    FiniteElementSpace fespace_nd(&mesh, &nd_fec);
    GridFunction dT(&fespace_nd);
    {
        DiscreteLinearOperator grad(&fespace, &fespace_nd);
        grad.AddDomainInterpolator(new GradientInterpolator);
        grad.Assemble();
        grad.Finalize();
        grad.Mult(x, dT);
        ofstream dT_ofs("dT.gf");
        dT_ofs.precision(8);
        dT.Save(dT_ofs);
    }

    {
        // 15. Save data in the ParaView format
        ParaViewDataCollection paraview_dc("ex0_potential", &mesh);
        paraview_dc.SetPrefixPath("ParaView");
        paraview_dc.SetLevelsOfDetail(order);
        paraview_dc.SetCycle(0);
        paraview_dc.SetDataFormat(VTKFormat::BINARY);
        paraview_dc.SetHighOrderOutput(true);
        paraview_dc.SetTime(0.0); // set the time
        paraview_dc.RegisterField("potential", &x);
        paraview_dc.RegisterField("sigma", &sigma_gf);
        paraview_dc.RegisterField("dT", &dT);
        paraview_dc.Save();
    }

   return 0;
}
