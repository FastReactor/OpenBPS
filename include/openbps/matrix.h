#ifndef MATRIX_H
#define MATRIX_H
#include <vector>
#include <memory>
#include <iostream>
#include <iomanip>
#include <initializer_list>
#include "xtensor/xarray.hpp"
#include "xtensor/xadapt.hpp"
#include "uncertainty.h"
#include "nuclide.h"
#include "chain.h"
#include "materials.h"
#include "reactions.h"

namespace openbps {
template <typename T>
class BaseMatrix;
class BaseMatrixIter;
class DecayMatrix;
class IterMatrix;
class CramMatrix;
//==============================================================================
// BaseMatrix description and implementation
//==============================================================================
template <typename T>
class BaseMatrix
{
public:
    //--------------------------------------------------------------------------
    //! Constructors, destructors, factory functions
    //! Initialize all matrix element to 0
    explicit BaseMatrix(size_t cols, size_t rows) : ncols_(cols), nrows_(rows)
    {allocate_memory();}
    //! Initialization by list
    BaseMatrix
    (std::initializer_list<std::initializer_list<T>> input) {
        size_t maxcols {1};
        this->nrows_ = input.size();
        for (auto it = input.begin(); it != input.end() &&
             it->size() > maxcols; it++)
            maxcols = it->size();
        ncols_ = maxcols;
        allocate_memory();
        size_t i{0}, j{0};
        for (auto v : input) {
            for (auto u: v) {
                data_[i][j] = u;
                j++;
            }
            i++;
            j = 0;
        }
    }
    //! Copy matrix
    BaseMatrix(const BaseMatrix& other) : nrows_(other.nrows_),
        ncols_ (other.ncols_) {
        allocate_memory();
        for (size_t i = 0; i < nrows_; i++) {
            std::copy(other.data_[i], other.data_[i] + ncols_,
                      this->data_[i]);
        }
    }

    //! '=' copy operator
    BaseMatrix& operator=(const BaseMatrix& other) {
        if (this != &other) {
            deallocate_memory();
            nrows_ = other.nrows_;
            ncols_ = other.ncols_;
            allocate_memory();
            for (size_t i = 0; i < nrows_; i++)
                std::copy(other.data_[i], other.data_[i] + ncols_,
                          this->data_[i]);
        }
        return *this;
    }
    //! Move matrix
    BaseMatrix(BaseMatrix &&other) : nrows_(other.nrows_),
        ncols_(other.ncols_) {
        allocate_memory();
        for (size_t i = 0; i < nrows_; i++) {
            data_[i] = other.data_[i];
        }
        other.data_ = nullptr;
        other.ncols_ = 0;
        other.nrows_ = 0;

    }

    //! '=' move operator
    BaseMatrix& operator=(BaseMatrix&& other) {
        if (this != &other) {
            deallocate_memory();
            this->nrows_ = other.nrows_;
            this->ncols_ = other.ncols_;
            allocate_memory();
            for (size_t i = 0; i < nrows_; i++)
                data_[i] = other.data_[i];
            other.data_ = nullptr;
            other.ncols_ = 0;
            other.nrows_ = 0;

        }
        return *this;
    }
    ~BaseMatrix() {deallocate_memory();}
    //--------------------------------------------------------------------------
    //! Methods and overloading operators
    //! Get a columns number
    //!
    size_t Numcols() {return ncols_;}
    //! Get a rows number
    //!
    size_t Numrows() const {return nrows_;}
    //! Indexing operator
    //!
    //! The indexing operator [][]
    T* operator [] (int i) {
        if (i >= nrows_)
            return nullptr;
        return data_[i];
    }
    //! M1 | M2
    const BaseMatrix& operator |(const BaseMatrix& M2) {
        if (this->nrows_ == M2.nrows_) {
            T* copyrow = new T[this->ncols_];
            for (size_t i = 0; i < this->nrows_; i++) {
                std::copy(this->data_[i], this->data_[i] + this->ncols_,
                          copyrow);
                delete[] this->data_[i] ;
                this->data_[i] = new T[this->ncols_ + M2.ncols_];
                std::copy(copyrow, copyrow + this->ncols_,
                          this->data_[i]);
                for (size_t j = 0; j < M2.ncols_; j++) {
                    this->data_[i][j + this->ncols_] = M2.data_[i][j];
                }
            }
            this->ncols_ += M2.ncols_;
        }
        return *this;
    }
    //! M1 + M2
    //!
    friend const BaseMatrix operator+(
            const BaseMatrix& M1,  const BaseMatrix& M2) {
        BaseMatrix M3(M1.ncols_, M1.nrows_);
        if (M1.ncols_ == M2.ncols_ && M1.nrows_ == M2.ncols_) {
            for (size_t i = 0; i < M3.ncols_; i++)
                for (size_t j = 0; j < M3.nrows_; j++)
                    M3.data_[i][j] = M1.data_[i][j] + M2.data_[i][j];
        }
        return M3;
    }
    //! Print out the matrix
    //!
    friend std::ostream& operator<< (std::ostream& out,
                                     const BaseMatrix& m1) {
        for (size_t i = 0; i < m1.nrows_; i++) {
            for (size_t j = 0; j < m1.ncols_; j++)
                out << std::setw(4)
                    << m1.data_[i][j] << "   ";
            out << std::endl;
        }
        return out;
    }

    //==========================================================================
    //! An iterator over BaseMatrix elements
    //==========================================================================

    class BaseMatrixIter
    {
    public:
      int indx_;  //!< An index
      int ncol_;  //!< Number of columns
      int nrow_;  //!< Number of columns

      BaseMatrixIter(BaseMatrix &mat, int ncol, int nrow, int indx)
        : indx_(indx), ncol_(ncol), nrow_(nrow), mat_(mat)
      {}

      bool operator==(const BaseMatrixIter &rhs) {return (indx_ == rhs.indx_);}

      bool operator!=(const BaseMatrixIter &rhs) {return !(*this == rhs);}

      T& operator*() {
          int irow {indx_ / ncol_};
          int icol {indx_ % ncol_};
          return mat_[irow][icol];}

      BaseMatrixIter& operator++()
      {
        if (indx_ < ncol_ * nrow_) {
          ++indx_;
          return *this;
        } else {
          indx_ = ncol_ * nrow_;
          return *this;
        }
      }

    protected:
      BaseMatrix& mat_;
    }; // class BaseMatrixIter

    //! For iteration start point
    //!
    BaseMatrixIter begin()
    {return BaseMatrixIter(*this, ncols_, nrows_, 0);}
    //! For iteration finish point
    //!
    BaseMatrixIter end()
    {return BaseMatrixIter(*this, ncols_, nrows_, ncols_ * nrows_);}
protected:
    //--------------------------------------------------------------------------
    //! Attributes
    //!
    size_t ncols_; //!< Columns number
    size_t nrows_; //!< Rows number
    T** data_;//!< Inner matrix data storage
    //--------------------------------------------------------------------------
    //! Methods
    //!
    //! Memory allocation
    void allocate_memory()
    {
        data_  = new T*[nrows_];
        for (size_t i = 0; i < nrows_; i++) {
            data_[i] = new T[ncols_];
            for (size_t j = 0; j < ncols_; j++)
                data_[i][j] = 0;
        }

    }

    //! Memory release
    void deallocate_memory()
    {
        if (data_ != nullptr){
            for (size_t i = 0; i < nrows_; i++)
                delete [] data_[i];
            delete [] data_;
            data_ = nullptr;
            nrows_ = 0;
            ncols_ = 0;
        }
    }
}; // class BaseMatrix

//==============================================================================
// DecayMatrix description
//==============================================================================

class DecayMatrix : public BaseMatrix<double> {
public:
    //--------------------------------------------------------------------------
    //! Constructors, destructors, factory functions
    explicit DecayMatrix(size_t dim) : BaseMatrix<double>(dim, dim) {}
    //--------------------------------------------------------------------------
    //! Methods
    //! Form a real decay nuclide matrix
    //!
    //!\param[in] chain with decay information and fill the data storage
    void form_matrixreal(Chain& chain);
    //! Form a deviation decay nuclide matrix for unceratanties analysis
    //!
    //!\param[in] chain with decay information and fill the data storage
    void form_matrixdev(Chain& chain);

}; //class DecayMatrix

//==============================================================================
// IterativeMatrix description
//==============================================================================

class IterMatrix : public DecayMatrix {
public:
    //--------------------------------------------------------------------------
    //! Constructors, destructors, factory functions
    explicit IterMatrix(size_t dim) : DecayMatrix(dim) {}
    IterMatrix(DecayMatrix& externmatrix) : DecayMatrix(externmatrix.Numrows()) {
        for (size_t i = 0; i < externmatrix.Numrows(); i++)
             for (size_t j = 0; j < externmatrix.Numrows(); j++)
                      this->data_[i][j] = externmatrix[i][j];

    }
    //--------------------------------------------------------------------------
    //! Methods
    //! Form a real decay nuclide matrix
    //!
    //!\param[in] chain with decay information and fill the data storage
    //!\param[in] material with nuclear concentration
    //!\return result a matrix with all transition for the material
    xt::xarray<double> matrixreal(Chain& chain,
                                  const Materials& mat);
    //! Form a deviation decay nuclide matrix for unceratanties analysis
    //!
    //!\param[in] chain with decay information and fill the data storage
    //!\param[in] material with nuclear concentration
    //!\return result a matrix with all transition for the material
    xt::xarray<double> matrixdev(Chain& chain,
                                 const Materials& mat);
    //! Form a real decay nuclide vector for all transition from nuclide
    //!
    //!\param[in] chain with decay information and fill the data storage
    //!\param[in] material with nuclear concentration
    //!\return result a vector with all transition for every nuclide:Real
    xt::xarray<double> sigp(Chain& chain,
                            const Materials& mat);
    //! Form a real decay nuclide vector for all transition from nuclide
    //!
    //!\param[in] chain with decay information and fill the data storage
    //!\param[in] material with nuclear concentration
    //!\return result a vector with all transition for every nuclide:Dev
    xt::xarray<double> dsigp(Chain& chain,
                             const Materials& mat);

}; //class IterMatrix

//==============================================================================
// ChebyshevMatrix description
//==============================================================================

class CramMatrix : public DecayMatrix {
public:
    //--------------------------------------------------------------------------
    //! Constructors, destructors, factory functions
    explicit CramMatrix(size_t dim) : DecayMatrix(dim) {}
    CramMatrix(DecayMatrix& externmatrix) : DecayMatrix(externmatrix.Numrows()) {
        for (size_t i = 0; i < externmatrix.Numrows(); i++)
             for (size_t j = 0; j < externmatrix.Numrows(); j++)
                      this->data_[i][j] = externmatrix[i][j];

    }
    //--------------------------------------------------------------------------
    //! Methods
    //! Form a real decay nuclide matrix
    //!
    //!\param[in] chain with decay information and fill the data storage
    //!\param[in] material with nuclear concentration
    //!\return result a matrix with all transition for the material
    xt::xarray<double> matrixreal(Chain& chain,
                                  const Materials& mat);
}; //class ChebyshevMatrix

//==============================================================================
// Non class methods
//==============================================================================

//! Get a concentration vector according to chain nuclide list
//!
//!\param[in] chainer chain instance
//!\param[in] nameconc nuclide names in current material
//!\param[in] ro vector with materials nuclide concentration
//!\param[in] isDev (by default false) wheter returns value contains Real or Dev
//!\return result nuclide concentration vector
xt::xarray<double> make_concentration(Chain& chainer, 
                                      const std::vector<std::string>& nameconc,
                                      const std::vector<udouble>& ro,
                                      bool isDev = false);

//! Find out power normalization coefficient to reaction-rate
//!
//!\param[in] mat Material to perform normalization
void power_normalization(Materials& mat);

} // namespace openbps
#endif // MATRIX_H
