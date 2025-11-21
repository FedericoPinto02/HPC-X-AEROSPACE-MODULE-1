#ifndef NSBSOLVER_INITIALIZER_HPP
#define NSBSOLVER_INITIALIZER_HPP

#include <memory>
#include <string>

#include "io/inputReader.hpp"
#include "core/Mesh.hpp"
#include "core/Fields.hpp"
#include "io/MuParserXAdapter.h"
#include "simulation/SimulationContext.hpp"

using TemporalFunc = std::function<double(double, double, double, double)>;
using SpatialFunc = std::function<double(double, double, double)>;

/**
 * @brief Class responsible for initializing the simulation context.
 */
class Initializer {
public:
    Initializer() = default;

    /**
     * @brief Setup and initialize the SimulationData from input data.
     */
    static SimulationData setup(const InputData& inputData);

    // // ---------------------------------------------------------------------
    // // Inizializzazione Campi Scalari (Field)
    // // ---------------------------------------------------------------------

    // /**
    //  * @brief Inizializza un Field scalare (e.g., k) da una funzione spaziale f(x, y, z).
    //  */
    // static Field initializeFieldFromSpatialFunc(
    //     const std::shared_ptr<Grid>& grid, 
    //     const SpatialFunc& func
    // );

    // /**
    //  * @brief Inizializza un Field scalare (e.g., pi) da una funzione temporale f(t, x, y, z) a un dato 'time'.
    //  */
    // static Field initializeFieldFromTemporalFunc(
    //     const double time,
    //     const std::shared_ptr<Grid>& grid, 
    //     const TemporalFunc& func
    // );

    // // ---------------------------------------------------------------------
    // // Inizializzazione Campi Vettoriali (VectorField)
    // // ---------------------------------------------------------------------

    // /**
    //  * @brief Inizializza un VectorField (e.g., campo iniziale u) da tre funzioni spaziali f(x, y, z).
    //  */
    // static VectorField initializeVectorFieldFromSpatialFunc(
    //     const std::shared_ptr<Grid>& grid,
    //     const SpatialFunc& func_u,
    //     const SpatialFunc& func_v,
    //     const SpatialFunc& func_w
    // );

    // /**
    //  * @brief Inizializza un VectorField (e.g., forze f) da tre funzioni temporali f(t, x, y, z) a un dato 'time'.
    //  */
    // static VectorField initializeVectorFieldFromTemporalFunc(
    //     const double time,
    //     const std::shared_ptr<Grid>& grid, // Corretto per coerenza
    //     const TemporalFunc& func_u,
    //     const TemporalFunc& func_v,
    //     const TemporalFunc& func_w
    // );

    // // ---------------------------------------------------------------------
    // // Aggiornamento Campi
    // // ---------------------------------------------------------------------

    // /**
    //  * @brief Aggiorna un VectorField esistente in-place (e.g., f, uBoundNew) al nuovo 'time'.
    //  */
    // static void updateVectorFieldWithTemporalFunc(
    //     const double time,
    //     VectorField& vec,
    //     const TemporalFunc& func_u,
    //     const TemporalFunc& func_v,
    //     const TemporalFunc& func_w
    // );

    static SpatialFunc makeSpatialFunc(const std::string& expr);

    static TemporalFunc makeTemporalFunc(const std::string& expr);

    static Field initializeFieldFromSpatialFunc(
        const std::shared_ptr<const Grid>& grid, // Aggiunto const
        const SpatialFunc& func
    );

    static Field initializeFieldFromTemporalFunc(
        const double time,
        const std::shared_ptr<const Grid>& grid, // Aggiunto const
        const TemporalFunc& func
    );

    static VectorField initializeVectorFieldFromSpatialFunc(
        const std::shared_ptr<const Grid>& grid, // Aggiunto const
        const SpatialFunc& func_u,
        const SpatialFunc& func_v,
        const SpatialFunc& func_w
    );

    static VectorField initializeVectorFieldFromTemporalFunc(
        const double time,
        const std::shared_ptr<const Grid>& grid, // Aggiunto const
        const TemporalFunc& func_u,
        const TemporalFunc& func_v,
        const TemporalFunc& func_w
    );

    static void updateVectorFieldWithTemporalFunc(
        const double time,
        VectorField& vec,
        const TemporalFunc& func_u,
        const TemporalFunc& func_v,
        const TemporalFunc& func_w
    );

    static Field initializeFieldFromExpr(const std::shared_ptr<const Grid>& grid, const std::string& expr);

    static VectorField initializeVectorFieldFromExpr(
        const std::shared_ptr<const Grid>& grid,
        const std::string& expr_u,
        const std::string& expr_v,
        const std::string& expr_w
    );

};

#endif // NSBSOLVER_INITIALIZER_HPP