#pragma once

#include <boost/qvm/mat_traits.hpp>
#include <boost/qvm/mat_operations.hpp>
#include <boost/qvm/map_mat_mat.hpp>
#include <ostream>

// code found at https://www.boost.org/doc/libs/1_73_0/libs/qvm/doc/html/index.html#mat_traits
namespace boost { namespace qvm {

  template <class T,int Rows,int Cols>
  struct mat {

    T a[Rows][Cols];

    template <class R>
    operator R() const {
      R r;
      assign(r,*this);
      return r;
    }

    void assign(std::vector<std::vector<T>>& array) {
        for (int i = 0; i < Rows; ++i) {
            for (int j = 0; j < Cols; ++j) {
                a[i][j] = array[i][j];
            }
        }
    }

    void assign_diagonal(std::vector<T>& array) {

        if (Rows != Cols)
            throw std::runtime_error("matrix must be square");

        for (int i = 0; i < Rows; ++i) {
            a[i][i] = array[i];
        }
    }

    std::string to_string() {
        std::string string = "";
        for (int i = 0; i < Rows; ++i) {
            std::string line = "";
            for (int j = 0; j < Cols; ++j) {
                line += std::to_string(a[i][j]) + " ";
            }
            line += "\n";
            string += line;
        }
        return string;
    }


  };

  template <class Matrix>
  struct mat_traits;

  template <class T,int Rows,int Cols>
  struct mat_traits< mat<T,Rows,Cols> > {

    typedef T scalar_type;
    static int const rows=Rows;
    static int const cols=Cols;

    template <int Row,int Col>
    static scalar_type read_element( mat<T,Rows,Cols> const & x ) {
      return x.a[Row][Col];
    }

    template <int Row,int Col>
    static scalar_type & write_element( mat<T,Rows,Cols> & x ) {
      return x.a[Row][Col];
    }

    static scalar_type read_element_idx( int row, int col, mat<T,Rows,Cols> const & x ) {
      return x.a[row][col];
    }

    static scalar_type & write_element_idx( int row, int col, mat<T,Rows,Cols> & x ) {
      return x.a[row][col];
    }

  };

} }