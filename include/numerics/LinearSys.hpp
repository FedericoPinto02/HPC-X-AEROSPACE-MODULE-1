#include <vector>
#include <core/TridiagMat.hpp>
#include <core/Fields.hpp>
#include <stdexcept>

enum class BoundaryType {Normal, Tangent};

class LinearSys {

private: 
    TridiagMat matA;
    std::vector<double> rhsC;
    std::vector<double> unknownX;
    BoundaryType boundaryType;

public:
     
    /**
     * @brief Constructor.
     * @param n Size of the linear system (n x n)
     * @param boundaryType Boundary condition type, Normal or Tangent (set Normal for pressure)
     */
    LinearSys(int n, BoundaryType boundaryType);


    /**
     * @brief Fill the linear system for Pressure variables
     */
    void fillSystemPressure(std::vector<double>& rhsIncomplete, Field& phi, const Axis direction);

    /**
     * @brief Fill the linear system for Velocity variables
     */
    void fillSystemVelocity(Field& porosity, std::vector<double>& rhsIncomplete, 
        VectorField& etaNew, VectorField& eta, VectorField& xi,
        const Axis direction, const size_t iStart, const size_t jStart, const size_t kStart);
    
    void ThomaSolver();
    
    /**
     * @brief Solve the linear system
     */
    void solveSystem();

   




};