#include <complex>
#include <iostream>

#include <mfem.hpp>

typedef std::complex<double> Complex;

const double PI = 3.14159265358979323846;
const double MU = 4 * PI * 1E-7;

const Complex II = Complex(0.0, 1.0);

const int XY_MODE = 0;
const int YX_MODE = 1;

const int Real = 0;
const int Imag = 1;

// Plane wave function for half-space model
class PlaneWaveFunction {
  public:
    void initialize(double freq, double z_top, double z_interface, double sig_air,
                    double sig_earth) {
        omega_ = 2 * PI * freq;

        z_tops_[0] = z_top;
        sig_[0] = sig_air;
        k_[0] = std::sqrt(-II * omega_ * MU * sig_[0]);

        z_tops_[1] = z_interface;
        sig_[1] = sig_earth;
        k_[1] = std::sqrt(-II * omega_ * MU * sig_[1]);

        double h0 = z_tops_[1] - z_tops_[0];
        Complex exp_k0h0 = std::exp(-k_[0] * h0);

        // Reflection coefficient
        Complex R = (k_[0] - k_[1]) / (k_[0] + k_[1]);
        Complex Rp1 = R * exp_k0h0;

        // Coefficients for Air (Layer 0)
        Complex b0 = 1.0 / (1.0 + Rp1 * exp_k0h0);
        Complex d0 = 1.0 / (1.0 - Rp1 * exp_k0h0);

        abcd_[0][0] = b0 * Rp1;
        abcd_[0][1] = b0;
        abcd_[0][2] = -d0 * Rp1;
        abcd_[0][3] = d0;

        // Coefficients for Earth (Layer 1)
        abcd_[1][0] = 0.0;
        abcd_[1][1] = abcd_[0][0] + abcd_[0][1] * exp_k0h0;
        abcd_[1][2] = 0.0;
        abcd_[1][3] = abcd_[0][2] + abcd_[0][3] * exp_k0h0;
    }

    void calculate_ehfield(double z, Complex eh[2]) const {
        int ilay = (z < z_tops_[1]) ? 0 : 1;

        Complex k = k_[ilay];
        Complex expm = std::exp(-k * (z - z_tops_[ilay]));
        Complex expp = (ilay == 0) ? std::exp(k * (z - z_tops_[1])) : Complex(0.0);

        Complex expmdz = -k * expm;
        Complex exppdz = k * expp;

        eh[0] = abcd_[ilay][0] * expp + abcd_[ilay][1] * expm;
        eh[1] = (abcd_[ilay][2] * exppdz + abcd_[ilay][3] * expmdz) / sig_[ilay];
    }

    void vector_value(int mode, int cc, const mfem::Vector &p, mfem::Vector &v) const {
        Complex eh[2];

        v.SetSize(3);
        v[0] = v[1] = v[2] = 0.0;

        calculate_ehfield(p[2], eh);
        if (cc == Real) {
            v[mode] = std::real(eh[mode]);
        } else {
            v[mode] = -std::imag(eh[mode]);
        }
    }

  private:
    double omega_, z_tops_[2], sig_[2];
    Complex k_[2], abcd_[2][4];
};

void forward() {
    // Create a simple half-space model with two cells
    double cell_size = 50000.0;
    mfem::Mesh serial_mesh = mfem::Mesh::MakeCartesian3D(1, 1, 2, mfem::Element::HEXAHEDRON,
                                                         cell_size, cell_size, 2 * cell_size);

    // Set mesh order for curved elements if needed, then call Mesh::Transform to
    // apply the true geometry
    int mesh_order = 1;
    serial_mesh.SetCurvature(mesh_order);

    // Create the parallel mesh
    mfem::ParMesh mesh(MPI_COMM_WORLD, serial_mesh);

    int myid = mesh.GetMyRank();

    if (myid == 0) {
        std::cout << "Number of Cells: " << mesh.GetGlobalNE() << std::endl;
    }

    // Polynomial degree
    int degree = 2;

    // Define finite element space for ND elements
    mfem::ND_FECollection fec(degree, 3, mfem::BasisType::GaussLobatto,
                              mfem::BasisType::IntegratedGLL);

    // Create the parallel finite element space
    mfem::ParFiniteElementSpace fespace(&mesh, &fec);

    if (myid == 0) {
        std::cout << "Number of DoFs: " << fespace.GlobalTrueVSize() << std::endl;
    }

    // Determine the list of essential boundary dofs
    mfem::Array<int> ess_bdr, ess_tdof_list;
    ess_bdr.SetSize(mesh.bdr_attributes.Max());
    ess_bdr = 1;
    fespace.GetEssentialTrueDofs(ess_bdr, ess_tdof_list);

    double freq = 0.01;
    double omega = 2 * PI * freq;

    // Define coefficients
    mfem::ConstantCoefficient mu_inv(1.0 / MU);
    mfem::FunctionCoefficient sigma([cell_size](const mfem::Vector &x) {
        if (x[2] <= cell_size) {
            return 1E-8; // Air
        } else {
            return 1E-2; // Half-space
        }
    });
    mfem::ProductCoefficient omega_sigma(omega, sigma);

    // Define the bilinear form for complex Maxwell's equations
    mfem::ParSesquilinearForm blf_a(&fespace, mfem::ComplexOperator::BLOCK_SYMMETRIC);
    blf_a.SetAssemblyLevel(mfem::AssemblyLevel::PARTIAL);
    blf_a.AddDomainIntegrator(new mfem::CurlCurlIntegrator(mu_inv), NULL);
    blf_a.AddDomainIntegrator(NULL, new mfem::VectorFEMassIntegrator(omega_sigma));
    blf_a.Assemble();

    // Define the plane wave boundary conditions
    PlaneWaveFunction pwf;
    pwf.initialize(freq, 0.0, cell_size, 1E-8, 1E-2);

    // Boundary function for real component of XY polarization
    mfem::VectorFunctionCoefficient e_real_xy(3, std::bind(&PlaneWaveFunction::vector_value, &pwf,
                                                           XY_MODE, Real, std::placeholders::_1,
                                                           std::placeholders::_2));

    // Boundary function for imaginary component of XY polarization
    mfem::VectorFunctionCoefficient e_imag_xy(3, std::bind(&PlaneWaveFunction::vector_value, &pwf,
                                                           XY_MODE, Imag, std::placeholders::_1,
                                                           std::placeholders::_2));

    // Boundary function for real component of YX polarization
    mfem::VectorFunctionCoefficient e_real_yx(3, std::bind(&PlaneWaveFunction::vector_value, &pwf,
                                                           YX_MODE, Real, std::placeholders::_1,
                                                           std::placeholders::_2));

    // Boundary function for imaginary component of YX polarization
    mfem::VectorFunctionCoefficient e_imag_yx(3, std::bind(&PlaneWaveFunction::vector_value, &pwf,
                                                           YX_MODE, Imag, std::placeholders::_1,
                                                           std::placeholders::_2));

    // Project boundary conditions
    mfem::ParComplexGridFunction e_xy(&fespace), e_yx(&fespace);
    e_xy.ProjectBdrCoefficientTangent(e_real_xy, e_imag_xy, ess_bdr);
    e_yx.ProjectBdrCoefficientTangent(e_real_yx, e_imag_yx, ess_bdr);

    // Define the linear form (right-hand side)
    mfem::ParComplexLinearForm lf_b(&fespace, mfem::ComplexOperator::BLOCK_SYMMETRIC);
    lf_b = 0.0;

    // Assemble the linear system
    mfem::OperatorHandle A;
    mfem::Vector x_xy, x_yx, b_xy, b_yx;
    blf_a.FormLinearSystem(ess_tdof_list, e_xy, lf_b, A, x_xy, b_xy);
    blf_a.FormLinearSystem(ess_tdof_list, e_yx, lf_b, A, x_yx, b_yx);

    // Define and assemble the preconditioner matrix
    mfem::ParBilinearForm blf_pc(&fespace);
    blf_pc.SetAssemblyLevel(mfem::AssemblyLevel::PARTIAL);
    blf_pc.AddDomainIntegrator(new mfem::CurlCurlIntegrator(mu_inv));
    blf_pc.AddDomainIntegrator(new mfem::VectorFEMassIntegrator(omega_sigma));
    blf_pc.Assemble();

    mfem::OperatorHandle B;
    blf_pc.FormSystemMatrix(ess_tdof_list, B);

    // Define the inner LOR preconditioner
    mfem::LORSolver<mfem::HypreAMS> inner_pc(blf_pc, ess_tdof_list);

    mfem::CGSolver inner_solver(mesh.GetComm());
    inner_solver.SetOperator(*B);
    inner_solver.SetPreconditioner(inner_pc);
    inner_solver.SetRelTol(1E-2);
    inner_solver.SetMaxIter(100);
    inner_solver.SetPrintLevel(1);

    mfem::Array<int> bdp_offsets;
    bdp_offsets.SetSize(3);
    bdp_offsets[0] = 0;
    bdp_offsets[1] = A->Height() / 2;
    bdp_offsets[2] = A->Height() / 2;
    bdp_offsets.PartialSum();

    // Define the block diagonal preconditioner
    mfem::BlockDiagonalPreconditioner bdp(bdp_offsets);
    bdp.SetDiagonalBlock(0, &inner_solver);
    bdp.SetDiagonalBlock(1, &inner_solver);

    mfem::FGMRESSolver fgmres(mesh.GetComm());
    fgmres.SetPreconditioner(bdp);
    fgmres.SetOperator(*A.Ptr());
    fgmres.SetRelTol(1E-8);
    fgmres.SetMaxIter(100);
    fgmres.SetPrintLevel(1);

    // Solve for XY and YX polarizations
    fgmres.Mult(b_xy, x_xy);
    blf_a.RecoverFEMSolution(x_xy, lf_b, e_xy);

    fgmres.Mult(b_yx, x_yx);
    blf_a.RecoverFEMSolution(x_yx, lf_b, e_yx);

    // Save the solutions for post-processing
    e_xy.real().SaveAsOne("e_xy_re.gf");
    e_xy.imag().SaveAsOne("e_xy_im.gf");

    e_yx.real().SaveAsOne("e_yx_re.gf");
    e_yx.imag().SaveAsOne("e_yx_im.gf");
}

int main(int argc, char **argv) {
    mfem::Mpi::Init(&argc, &argv);
    mfem::Hypre::Init();

    forward();

    return 0;
}
