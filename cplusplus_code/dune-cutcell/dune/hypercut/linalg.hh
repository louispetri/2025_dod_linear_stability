#ifndef DUNE_HYPERCUT_LINALG_HH
#define DUNE_HYPERCUT_LINALG_HH

#include <dune/common/fmatrix.hh>

namespace Dune {

    template< class K >
    class VLAVectorView
    {
        using value_type = K;
        using size_type = std::size_t;

    private:
        value_type* begin_;
        value_type* end_;

    public:
        VLAVectorView(value_type* begin, value_type* end)
            : begin_(begin)
            , end_(end)
        {

        }

        //==== make this thing a vector
        size_type size() const
        {
            return std::distance(begin_, end_);
        }

        K & operator[](size_type i) {
            DUNE_ASSERT_BOUNDS(i < size());
            return begin_[i];
        }
        const K & operator[](size_type i) const {
            DUNE_ASSERT_BOUNDS(i < size());
            return begin_[i];
        }

        //! return pointer to underlying array
        K* data() noexcept
        {
            return begin_;
        }

        //! return pointer to underlying array
        const K* data() const noexcept
        {
            return begin_;
        }
    };

    template< class K >
    class ConstVLAVectorView
    {
    public:
        using value_type = K;
        using size_type = std::size_t;

    private:
        const value_type* begin_;
        const value_type* end_;

    public:
        ConstVLAVectorView(const value_type* begin, const value_type* end)
            : begin_(begin)
            , end_(end)
        {

        }

        //==== make this thing a vector
        size_type size() const
        {
            return std::distance(begin_, end_);
        }

        const K & operator[](size_type i) const {
            DUNE_ASSERT_BOUNDS(i < size());
            return begin_[i];
        }

        //! return pointer to underlying array
        const K* data() const noexcept
        {
            return begin_;
        }
    };

    template<class K>
    class VLAMatrix
    {
        std::vector< K > data_;
        std::size_t row_count;
        std::size_t column_count;

        typedef DenseMatrix< VLAMatrix<K> > Base;

    public:
        using value_type = K;
        using size_type = std::size_t;
        using row_type = VLAVectorView<K>;
        using const_row_type = ConstVLAVectorView<K>;

        //===== constructors
        //! \brief Default constructor
        VLAMatrix () :
            row_count(0),
            column_count(0)
        {}

        //! \brief Constructor initializing the whole matrix with a scalar
        VLAMatrix (size_type r, size_type c, value_type v = value_type() ) :
            data_(r * c, v ),
            row_count(r),
            column_count(c)
        {}

        void resize (size_type r, size_type c, value_type v = value_type() )
        {
            row_count = r;
            column_count = c;
            data_.resize(r * c, v);
        }

        //===== assignment
        // General assignment with resizing
        template <typename T, typename = std::enable_if_t<!Dune::IsNumber<T>::value>>
        VLAMatrix& operator=(T const& rhs) {
            data_.resize(rhs.N() * rhs.M());
            Base::operator=(rhs);
            return *this;
        }

        // Specialisation: scalar assignment (no resizing)
        template <typename T, typename = std::enable_if_t<Dune::IsNumber<T>::value>>
        VLAMatrix& operator=(T scalar) {
            std::fill(data_.begin(), data_.end(), scalar);
            return *this;
        }

        void transposed(VLAMatrix<K>& AT) const
        {
            AT.resize(column_count, row_count);

            for( size_type i = 0; i < mat_rows(); ++i ) {
                for( size_type j = 0; j < mat_cols(); ++j ) {
                    AT[j][i] = (*this)[i][j];
                }
            }
        }

        // make this thing a matrix
        size_type mat_rows() const
        {
            return row_count;
        }

        size_type mat_cols() const
        {
            return column_count;
        }

        row_type operator[](size_type i)
        {
            DUNE_ASSERT_BOUNDS(i < row_count);
            value_type* ptr = data_.data() + i * column_count;
            return row_type(ptr, ptr + column_count);
        }

        const_row_type operator[](size_type i) const
        {
            DUNE_ASSERT_BOUNDS(i < row_count);
            const value_type* ptr = data_.data() + i * column_count;
            return const_row_type(ptr, ptr + column_count);
        }
    };

    template<class K>
    std::ostream& operator<<(std::ostream& out, const VLAMatrix<K>& m)
    {
        for (std::size_t i = 0; i < m.mat_rows(); ++i) {
            for (std::size_t j = 0; j < m.mat_cols(); ++j) {
                out << m[i][j] << " ";
            }

            out << std::endl;
        }

        return out;
    }
}

namespace Dune::Hypercut {
    extern "C" void dgesvd_(const char* jobu, const char* jobvt,
                    const long int* m, const long int* n, double* a, const long int* lda,
                    double* s, double* u, const long int* ldu,
                    double* vt, const long int* ldvt,
                    double* work, const long int* lwork, long int* info);

    extern "C" void dscal_(const long int* n, double* da, double* dx, const long int* incx);

    extern "C" void dgemm_(const char* transa, const char* transb,
                           const long int* m, const long int* n, const long int* k,
                           double* alpha, double* a, const long int* lda,
                           double* b, const long int* ldb, double* beta,
                           double* c, const long int* ldc);

    template <int dim, typename K>
    void pseudoInverse(const Dune::FieldMatrix<K, dim, dim>& matrix, Dune::FieldMatrix<K, dim, dim>& result)
    {
        const long int N = dim ;
        const char jobu = 'A';
        const char jobvt = 'A';

        constexpr bool isKLapackType = std::is_same_v<K, double> || std::is_same_v<K, float>;
        using LapackNumType = std::conditional_t<isKLapackType, K, double>;

        LapackNumType matrixVector[dim * dim];

        int flatIndex = 0;

        for (int i = 0; i < dim; ++i) {
            for (int j = 0; j < dim; ++j, ++flatIndex) {
                matrixVector[flatIndex] = matrix[j][i];
            }
        }

        LapackNumType singularValues[dim];
        LapackNumType U[dim * dim];
        LapackNumType VT[dim * dim];

        LapackNumType work[10 * dim];
        const long int lwork = 10 * dim;

        long int info = 0;

        dgesvd_(&jobu, &jobvt, &N, &N, &matrixVector[0], &N, &singularValues[0], &U[0], &N, &VT[0], &N, &work[0], &lwork, &info);

        const long int incx = 1;

        for (int i = 0; i < N; ++i) {
            LapackNumType invSinVal = 0.0;

            if (singularValues[i] > 1e-13) {
                invSinVal = 1.0 / singularValues[i];
            }

            dscal_(&N, &invSinVal, &U[i * N], &incx);
        }

        LapackNumType psInv[dim * dim];

        const char transa = 'T';
        const char transb = 'T';

        double alpha = 1.0;
        double beta = 0.0;

        dgemm_(&transa, &transb, &N, &N, &N, &alpha, &VT[0],&N, &U[0], &N, &beta, &psInv[0], &N);

        flatIndex = 0;

        for (int i = 0; i < dim; ++i) {
            for (int j = 0; j < dim; ++j, ++flatIndex) {
                result[j][i] = psInv[flatIndex];
            }
        }
    }

    template <int dim, typename K>
    void rangeSpace(const Dune::FieldMatrix<K, dim, dim>& matrix, Dune::FieldMatrix<K, dim, dim>& result)
    {
        const long int N = dim ;
        const char jobu = 'A';
        const char jobvt = 'A';

        constexpr bool isKLapackType = std::is_same_v<K,double> || std::is_same_v<K,float>;
        using LapackNumType = std::conditional_t<isKLapackType, K, double>;

        LapackNumType matrixVector[dim * dim];

        int flatIndex = 0;

        for (int i = 0; i < dim; ++i) {
            for (int j = 0; j < dim; ++j, ++flatIndex) {
                matrixVector[flatIndex] = matrix[j][i];
            }
        }

        LapackNumType singularValues[dim];
        LapackNumType U[dim * dim];
        LapackNumType VT[dim * dim];

        LapackNumType work[10 * dim];
        const long int lwork = 10 * dim;

        long int info = 0;

        dgesvd_(&jobu, &jobvt, &N, &N, &matrixVector[0], &N, &singularValues[0], &U[0], &N, &VT[0], &N, &work[0], &lwork, &info);

        for (int i = 0; i < dim; ++i) {
            if (singularValues[i] < 1e-13) {
                break;
            }

            for (int j = 0; j < dim; ++j, ++flatIndex) {
                result[j][i] = U[flatIndex];
            }
        }
    }

    // assumes that the matrix is symmetric
    template<int dim, class K>
    void removeNegativeEigenvalues(Dune::FieldMatrix<K, dim, dim>& matrix, K epsilon = 1e-13)
    {
        using Matrix = Dune::FieldMatrix<K, dim, dim>;
        Matrix eigenvectors(0.0);
        Dune::FieldVector<K, dim> eigenvalues(0.0);

        Dune::FMatrixHelp::eigenValuesVectors(matrix, eigenvalues, eigenvectors);

        Matrix matrixNeg(0.0);

        for (std::size_t i = 0; i < matrixNeg.N(); ++i) {
            if (eigenvalues[i] < -epsilon) {
                matrixNeg[i][i] = eigenvalues[i];
            }
        }

        matrix = eigenvectors.transposed() * matrixNeg * eigenvectors;
    }
}

#endif