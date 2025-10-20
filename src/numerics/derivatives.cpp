#include "numerics/derivatives.hpp"
#include "core/Mesh.hpp"
#include <stdexcept>

void Derivatives::computeGradient(const Field &field, VectorField &gradient) {
    computeDx(field, gradient.x());
    computeDy(field, gradient.y());
    computeDz(field, gradient.z());
    // todo - may be optimized in terms of locality: tiling ??
}

void Derivatives::computeDx(const Field &field, Field &dx) {

    // safety checks
    // todo - are they really needed ?? wait profiling
    if (!field.getGrid()) {
        throw std::runtime_error("input field has null grid.");
    }
    if (!dx.getGrid()) {
        throw std::runtime_error("output field (dx) has null grid.");
    }
    if (field.getGrid() != dx.getGrid()) {
        throw std::runtime_error("field and dx must share the same grid.");
    }
    
    std::shared_ptr<const Grid> grid = field.getGrid();
    const size_t Nx = grid->Nx;
    const size_t Ny = grid->Ny;
    const size_t Nz = grid->Nz;
    const double mul = 0.5/grid->dx;

    if (Nx < 3) {
        throw std::runtime_error("grid size Nx must be at least 3 to compute central derivative.");
    }
    if (grid->dx <= 0.0) {
        throw std::runtime_error("grid spacing dx must be positive.");
    }
    if (Nx == 0 || Ny == 0 || Nz == 0) {
        throw std::runtime_error("grid dimensions must be > 0.");
    }

    for (size_t k=0 ; k<Nz ; k++)
        for (size_t j=0 ; j<Ny ; j++)
            for (size_t i=1 ; i<(Nx-1) ; i++)
                {
                    double dfdx = 0.0;
                    dfdx = (field(i + 1, j, k) - field(i - 1, j, k))*mul;
                    dx(i,j,k) = dfdx;
                }
}

void Derivatives::computeDy(const Field &field, Field &dy) {

    // safety checks
    // todo - are they really needed ?? wait profiling
    if (!field.getGrid()) {
        throw std::runtime_error("input field has null grid.");
    }
    if (!dy.getGrid()) {
        throw std::runtime_error("output field (dy) has null grid.");
    }
    if (field.getGrid() != dy.getGrid()) {
        throw std::runtime_error("field and dx must share the same grid.");
    }
    
    std::shared_ptr<const Grid> grid = field.getGrid();
    const size_t Nx = grid->Nx;
    const size_t Ny = grid->Ny;
    const size_t Nz = grid->Nz;
    const double mul = 0.5/grid->dy;

    if (Ny < 3) {
        throw std::runtime_error("grid size Ny must be at least 3 to compute central derivative.");
    }
    if (grid->dy <= 0.0) {
        throw std::runtime_error("grid spacing dx must be positive.");
    }
    if (Nx == 0 || Ny == 0 || Nz == 0) {
        throw std::runtime_error("grid dimensions must be > 0.");
    }

    for (size_t k=0 ; k<Nz ; k++)
        for (size_t i=0 ; i<Nx ; i++)
            for (size_t j=1 ; j<(Ny-1) ; j++)
                {
                    double dfdy = 0.0;
                    dfdy = (field(i, j + 1, k) - field(i, j - 1, k)) * mul;
                    dy(i,j,k) = dfdy;
                }
}

void Derivatives::computeDz(const Field &field, Field &dz) {
    
    // safety checks
    // todo - are they really needed ?? wait profiling
    if (!field.getGrid()) {
        throw std::runtime_error("input field has null grid.");
    }
    if (!dz.getGrid()) {
        throw std::runtime_error("output field (dz) has null grid.");
    }
    if (field.getGrid() != dz.getGrid()) {
        throw std::runtime_error("field and dz must share the same grid.");
    }
    
    std::shared_ptr<const Grid> grid = field.getGrid();
    const size_t Nx = grid->Nx;
    const size_t Ny = grid->Ny;
    const size_t Nz = grid->Nz;
    const double mul = 0.5/grid->dz;

    if (Nz < 3) {
        throw std::runtime_error("grid size Nz must be at least 3 to compute central derivative.");
    }
    if (grid->dz <= 0.0) {
        throw std::runtime_error("grid spacing dz must be positive.");
    }
    if (Nx == 0 || Ny == 0 || Nz == 0) {
        throw std::runtime_error("grid dimensions must be > 0.");
    }

    for (size_t i=0 ; i<Nx ; i++)
        for (size_t j=0 ; j<Ny ; j++)
            for (size_t k=1 ; k<(Nz-1) ; k++)
                {
                    double dfdz = 0.0;
                    dfdz = (field(i, j, k+1) - field(i, j, k-1)) * mul;
                    dz(i,j,k) = dfdz;
                }
}

void Derivatives::computeDivergence(const Field &field, Field &divergence) {
    Field tmp;
    tmp.setup(field.getGrid(),
                std::vector<Field::Scalar>(field.getGrid()->size(), 0.0));

    divergence.reset();
    computeDx(field, tmp);
    divergence.add(tmp);
    computeDy(field, tmp);
    divergence.add(tmp);
    computeDz(field, tmp);
    divergence.add(tmp);
}

void Derivatives::computeHessianDiag(const Field &field, VectorField &hessianDiag) {
    computeDx(field, hessianDiag.x());
    computeDy(field, hessianDiag.y());
    computeDz(field, hessianDiag.z());
    // todo - may be optimized in terms of locality: tiling ??
}

void Derivatives::computeDxx(const Field &field, Field &dxx) {

    // safety checks
    
    // todo - are they really needed ?? wait profiling
    if (!field.getGrid()) {
        throw std::runtime_error("input field has null grid.");
    }
    if (!dxx.getGrid()) {
        throw std::runtime_error("output field (dx) has null grid.");
    }
    if (field.getGrid() != dxx.getGrid()) {
        throw std::runtime_error("field and dx must share the same grid.");
    }
    
    std::shared_ptr<const Grid> grid = field.getGrid();
    const size_t Nx = grid->Nx;
    const size_t Ny = grid->Ny;
    const size_t Nz = grid->Nz;
    const double mul = 1.0 / (grid->dx * grid->dx);

    if (Nx < 3) {
        throw std::runtime_error("grid size Nx must be at least 3 to compute second derivative.");
    }
    if (grid->dx <= 0.0) {
        throw std::runtime_error("grid spacing dx must be positive.");
    }
    if (Nx == 0 || Ny == 0 || Nz == 0) {
        throw std::runtime_error("grid dimensions must be > 0.");
    }

    for (size_t k=0 ; k<Nz ; k++)
        for (size_t j=0 ; j<Ny ; j++)
            for (size_t i=1 ; i<(Nx-1) ; i++)
                {
                    double dfdxx = 0.0;
                    dfdxx = (field(i + 1, j, k) + field(i - 1, j, k) - 2.0 * field(i, j, k))*mul;
                    dxx(i,j,k) = dfdxx;
                }
}

void Derivatives::computeDyy(const Field &field, Field &dyy) {

    // safety checks
    // todo - are they really needed ?? wait profiling
    if (!field.getGrid()) {
        throw std::runtime_error("input field has null grid.");
    }
    if (!dyy.getGrid()) {
        throw std::runtime_error("output field (dy) has null grid.");
    }
    if (field.getGrid() != dyy.getGrid()) {
        throw std::runtime_error("field and dy must share the same grid.");
    }
    
    std::shared_ptr<const Grid> grid = field.getGrid();
    const size_t Nx = grid->Nx;
    const size_t Ny = grid->Ny;
    const size_t Nz = grid->Nz;
    const double mul = 1.0 / (grid->dy * grid->dy);

    if (Ny < 3) {
        throw std::runtime_error("grid size Ny must be at least 3 to compute second derivative.");
    }
    if (grid->dy <= 0.0) {
        throw std::runtime_error("grid spacing dy must be positive.");
    }
    if (Nx == 0 || Ny == 0 || Nz == 0) {
        throw std::runtime_error("grid dimensions must be > 0.");
    }

    for (size_t k=0 ; k<Nz ; k++)
        for (size_t j=1 ; j<(Ny-1) ; j++)
            for (size_t i=0 ; i<Nx ; i++)
                {
                    double dfdyy = 0.0;
                    dfdyy = (field(i, j+1, k) + field(i, j-1, k) - 2.0 * field(i, j, k))*mul;
                    dyy(i,j,k) = dfdyy;
                }
}

void Derivatives::computeDzz(const Field &field, Field &dzz) {

    // safety checks
    // todo - are they really needed ?? wait profiling
    if (!field.getGrid()) {
        throw std::runtime_error("input field has null grid.");
    }
    if (!dzz.getGrid()) {
        throw std::runtime_error("output field (dz) has null grid.");
    }
    if (field.getGrid() != dzz.getGrid()) {
        throw std::runtime_error("field and dz must share the same grid.");
    }
    
    std::shared_ptr<const Grid> grid = field.getGrid();
    const size_t Nx = grid->Nx;
    const size_t Ny = grid->Ny;
    const size_t Nz = grid->Nz;
    const double mul = 1.0 / (grid->dz * grid->dz);

    if (Nz < 3) {
        throw std::runtime_error("grid size Nz must be at least 3 to compute second derivative.");
    }
    if (grid->dz <= 0.0) {
        throw std::runtime_error("grid spacing dz must be positive.");
    }
    if (Nx == 0 || Ny == 0 || Nz == 0) {
        throw std::runtime_error("grid dimensions must be > 0.");
    }

    for (size_t k=1 ; k<(Nz-1) ; k++)
        for (size_t j=0 ; j<Ny ; j++)
            for (size_t i=0 ; i<Nx ; i++)
                {
                    double dfdzz = 0.0;
                    dfdzz = (field(i, j, k+1) + field(i, j, k-1) - 2.0 * field(i, j, k))*mul;
                    dzz(i,j,k) = dfdzz;
                }
}