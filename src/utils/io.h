#ifndef IO_H
#define IO_H

#include <cstddef>
#include <string>
#include <vector>
#include "Curve.h"
#include "geom.h"

namespace io {
    /**
     * \brief Read the curve coordinates from a file, skipping the header.
     * \param filename The name of the file to read.
     * \param curve The resulting curve.
     * \param header_size The number of lines of the file header.
     */
    void read_curve(std::string const& filename, Curve& curve,
        std::size_t header_size = 0);

    /**
     * \brief Read the curve coordinates from a file, skipping the header, and
     * return the resulting curve.
     * \param filename The name of the file to read.
     * \param header_size The number of lines of the file header.
     * \return The curve.
     */
    Curve read_curve(std::string const& filename, std::size_t header_size = 0);

    /**
     * \brief Read all the curves from dataset.txt in a directory.
     * \param directory The data directory.
     * \param header_size The number of lines of the file header.
     * \return The curves.
     */
    std::vector<Curve> read_curves(std::string const& directory,
        std::size_t header_size = 0);

    /**
     * \brief Export the data in a unified format into a file.
     * \param filename The output file.
     * \param points The curve / point sequence to write.
     */
    void export_points(std::string const& filename, Points const& points);

    /**
     * @brief Writes a warping path to file
     * 
     * @param filename 
     * @param edges 
     * @param n 
     * @param m 
     */
    void write_path(const std::string& filename, std::vector<std::vector<double>>& edges, int n, int m);
}
#endif
