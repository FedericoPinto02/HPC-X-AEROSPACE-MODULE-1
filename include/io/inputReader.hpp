#ifndef INPUTREADER_HPP
#define INPUTREADER_HPP

#include <iostream>
#include <string>
#include <extern/json.hpp>



struct InputData {
    struct MeshData {
        int nx, ny, nz;
        double dx, dy, dz;
    } mesh;

    struct PhysicsData {
        double viscosity;
    } physics;

    struct InitialConditions {
        double u0, v0, w0, p0;
    } initialConditions;

    struct TimeData {
        double dt;
        double t_end;
    } time;
};
class InputReader {

    public:
        InputReader() {};

        InputData readAndSetInput(const std::string& filename) {}

    private:
        std::string filename;


};

#endif // INPUTREADER_HPP