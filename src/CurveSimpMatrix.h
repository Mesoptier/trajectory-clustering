#pragma once

#include "Curve.h"
#include "IntegralFrechet/IntegralFrechet.h"
#include <fstream>

using Curves = std::vector<Curve>;

class CurveSimpMatrix {

    private:
        std::vector<std::vector<distance_t>> matrix;

    public:
        CurveSimpMatrix(std::string file_name) {
            std::fstream file;
            file.open(file_name);
            read(file);
        }

        CurveSimpMatrix(Curves const& curves, Curves const& simplifications, bool speed_up) {
            using szt = Curves::size_type;
            
            if (curves.size() != simplifications.size())
                throw std::runtime_error("There must be the same number of curves and simplifications in CurveSimpMatrix constructor");

            matrix = std::vector<std::vector<distance_t>>();

            for (szt i = 0; i < curves.size(); ++i) {
                matrix.push_back(std::vector<distance_t>());
                for (szt j = 0; j < curves.size(); ++j) {
                    if (!speed_up) {
                            matrix.back().push_back(
                            IntegralFrechet(curves[i], simplifications[j], ParamMetric::LInfinity_NoShortcuts, 1, nullptr)
                            .compute_matching()
                            .cost
                        );
                    } else {

                        const auto result_alt = IntegralFrechet(curves[i].simplify(true), simplifications[j], ParamMetric::LInfinity_NoShortcuts, 10, nullptr).compute_matching();
                        const auto band = MatchingBand(curves[i], simplifications[j], result_alt.matching, 1);
                        const auto result = IntegralFrechet(curves[i], simplifications[j], ParamMetric::LInfinity_NoShortcuts, 1, &band).compute_matching();

                        matrix.back().push_back(result.cost);

                    }
                }
            }

        }

        distance_t at(std::size_t i, std::size_t j) const {
            return matrix[i][j];
        }

        void write(std::ofstream& file) {
            for (auto row: matrix) {
                for (auto distance: row) {
                    file << distance << " ";
                }
                file << "\n";
            }
        }

        void read(std::fstream& file) {
            matrix = std::vector<std::vector<distance_t>>();

            for (std::string line; std::getline(file, line);) {
                matrix.push_back(std::vector<distance_t>());
                std::istringstream tokenStream(line);
                std::string distance;
                while (std::getline(tokenStream, distance, ' ')) {
                    matrix.back().push_back(std::stod(distance));
                }
            }
        }
};
