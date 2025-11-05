#ifndef NBSOLVER_MUPARSERXADAPTER_H
#define NBSOLVER_MUPARSERXADAPTER_H

#include <cmath>
#include <functional>
#include <iostream>
#include <memory>
#include <string>

#include "mpParser.h"

namespace MuParserXAdapter {

    /**
     * @brief Creates a std::function with x, y, z inputs from a muparserx expression.
     * @param expr the muparserx expression as a string (e.g., "sin(pi*x) * cos(pi*y) * exp(-z)")
     * @return the C++ function that can be called with (x,y,z) to evaluate the expression;
     * zero function if parsing fails
     */
    inline std::function<double(double x, double y, double z)>
    createFunction(const std::string &expr) {
        /*
         * Bind parser and (x,y,z) variables to the returned lambda via shared pointers,
         * so that they remain valid even after this procedure returns.
         * */
        auto parser = std::make_shared<mup::ParserX>();
        auto x_var = std::make_shared<mup::Value>(0.0);
        auto y_var = std::make_shared<mup::Value>(0.0);
        auto z_var = std::make_shared<mup::Value>(0.0);

        try {
            // Define and bind the (x,y,z) expression variables to actual C++ variables
            parser->DefineVar("x", mup::Variable(x_var.get()));
            parser->DefineVar("y", mup::Variable(y_var.get()));
            parser->DefineVar("z", mup::Variable(z_var.get()));

            parser->SetExpr(expr);

            // Create the lambda "adapter"
            auto func = [parser, x_var, y_var, z_var](double x, double y, double z) -> double {
                // Set muparser x,y,z variables with input values
                *x_var = x;
                *y_var = y;
                *z_var = z;
                return parser->Eval().GetFloat();
            };
            return func;
        } catch (mup::ParserError &e) {
            std::cerr << "muparserx error: " << e.GetMsg() << std::endl;
            std::cerr << "Using zero function as fallback." << std::endl;
            return [](double, double, double) { return 0.0; };
        }
    }


    /**
     * @brief Creates a std::function with t, x, y, z inputs from a muparserx expression.
     * @param expr the muparserx expression as a string (e.g., "t * sin(pi*x) * cos(pi*y) * exp(-z)")
     * @return the C++ function that can be called with (t,x,y,z) to evaluate the expression;
     * zero function if parsing fails
     */
    inline std::function<double(double t, double x, double y, double z)>
    createTimeFunction(const std::string &expr) {
        /*
         * Bind parser and (t,x,y,z) variables to the returned lambda via shared pointers,
         * so that they remain valid even after this procedure returns.
         * */
        auto parser = std::make_shared<mup::ParserX>();
        auto t_var = std::make_shared<mup::Value>(0.0);
        auto x_var = std::make_shared<mup::Value>(0.0);
        auto y_var = std::make_shared<mup::Value>(0.0);
        auto z_var = std::make_shared<mup::Value>(0.0);

        try {
            // Define and bind the (t,x,y,z) expression variables to actual C++ variables
            parser->DefineVar("t", mup::Variable(t_var.get()));
            parser->DefineVar("x", mup::Variable(x_var.get()));
            parser->DefineVar("y", mup::Variable(y_var.get()));
            parser->DefineVar("z", mup::Variable(z_var.get()));

            parser->SetExpr(expr);

            // Create the lambda "adapter"
            auto func = [parser, t_var, x_var, y_var, z_var](double t, double x, double y, double z) -> double {
                // Set muparser t,x,y,z variables with input values
                *t_var = t;
                *x_var = x;
                *y_var = y;
                *z_var = z;
                return parser->Eval().GetFloat();
            };
            return func;
        } catch (mup::ParserError &e) {
            std::cerr << "muparserx error: " << e.GetMsg() << std::endl;
            std::cerr << "Using zero function as fallback." << std::endl;
            return [](double, double, double, double) { return 0.0; };
        }
    }

} // namespace MuParserXAdapter


#endif //NBSOLVER_MUPARSERXADAPTER_H
