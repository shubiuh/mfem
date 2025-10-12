#include "mfem.hpp" // Include the MFEM header

using namespace mfem;

int main(int argc, char* argv[]) {
    // Set up a simple mesh and finite element space
    int order = 2;
    int dim = 3;
    int n = 1;

    // Create a 3D mesh
    Mesh mesh = Mesh::MakeCartesian3D(n, n, n, Element::TETRAHEDRON, 1.0, 1.0, 1.0);

    // Create finite element collections and spaces
    H1_FECollection h1_fec(order, dim);
    FiniteElementSpace h1_fespace(&mesh, &h1_fec);

    // Create a GridFunction to hold the solution
    GridFunction h1_x(&h1_fespace);

    // Define a function coefficient
    FunctionCoefficient funcCoef([](const Vector &x) { return x[0] + x[1] + x[2]; });

    // Project the coefficient onto the GridFunction
    h1_x.ProjectCoefficient(funcCoef);

    // Create a GradientGridFunctionCoefficient
    GradientGridFunctionCoefficient h1_xCoef(&h1_x);

    // Test the evaluation of the gradient at a point
    Vector grad(dim);
    double coords[3] = {0.5, 0.5, 0.5}; // Example coordinates in the mesh
    IntegrationPoint ip = coords; // Create IntegrationPoint with coordinates
    h1_xCoef.Eval(grad, h1_xCoef, ip);

    // Output the gradient
    std::cout << "Gradient at point (" << ip.x << ", " << ip.y << ", " << ip.z << "): "
              << grad << std::endl;

    return 0;
}