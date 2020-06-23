#pragma once

#include "Curve.h"
#include "IntegralFrechet/IntegralFrechet.h"
#include <fstream>

using Curves = std::vector<Curve>;

class CurveSimpMatrix {

    private:
        std::vector<std::vector<distance_t>> matrix;

    public:
        CurveSimpMatrix() {};

        CurveSimpMatrix(std::string file_name) {
            std::fstream file;
            file.open(file_name);
            read(file);
            file.close();
        }

        CurveSimpMatrix(Curves const& curves, Curves const& simplifications, distance_t(*dist_func)(Curve, Curve)) {
            using szt = Curves::size_type;
            if (curves.size() != simplifications.size())
                throw std::runtime_error("There must be the same number of curves and simplifications in CurveSimpMatrix constructor");

            matrix = std::vector<std::vector<distance_t>>();


            for (szt i = 0; i < curves.size(); ++i) {
                matrix.push_back(std::vector<distance_t>());
                for (szt j = 0; j < curves.size(); ++j) {
                    if ((i * curves.size() + j) % 25 == 0)
                        std::cout << i * curves.size() + j << " / " << curves.size() * curves.size() << std::endl;
                    // if (i == j) {
                    //     matrix.back().push_back(0);
                    // }
                    // else {
                        matrix.back().push_back(dist_func(curves[i], simplifications[j]));
                    // }
                }
            }
        }

        distance_t at(std::size_t i, std::size_t j) const {
            return matrix[i][j];
        }

        void set(std::size_t i, std::size_t j, distance_t value) {
            matrix[i][j] = value;
        }

        void write(std::string path) {
            std::ofstream file;
            file.open(path, std::fstream::out | std::fstream::trunc);
            for (auto row: matrix) {
                for (auto distance: row) {
                    file << distance << " ";
                }
                file << "\n";
            }
            file.close();
        }

        void read(std::fstream& file) {
            matrix = std::vector<std::vector<distance_t>>();
            std::string line;
            int i = 0;
            while (std::getline(file, line)) {
                matrix.push_back(std::vector<distance_t>());
                std::istringstream tokenStream(line);
                std::string distance;
                while (std::getline(tokenStream, distance, ' ')) {
                    matrix.back().push_back(std::stod(distance));
                }
            }
        }
};
