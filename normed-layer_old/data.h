#ifndef __DATA__
#define __DATA__
#include <malloc.h>
#include <fstream>
#include <cmath>
#include <vector>

// Dataset
// T: element type
// A: alignment, default = 128 bits
#define MATRIX_ALIGN 1


template <typename T, int A = MATRIX_ALIGN>
class Dataset
{
    int dim;
    long long N;
    size_t stride;
    char *dims;
public:
    typedef T value_type;
    static const int ALIGN = A;
    void reset (int _dim, long long _N)
    {
        //BOOST_ASSERT((ALIGN % sizeof(T)) == 0);
        dim = _dim;
        N = _N;
        stride = dim * sizeof(T) + ALIGN - 1;
        stride = stride / ALIGN * ALIGN;
        if (dims != NULL) delete[] dims;
        dims = (char *)memalign(ALIGN, N * stride); // SSE instruction needs data to be aligned
        std::fill(dims, dims + N * stride, 0);
    }

    void free (void) {
        dim = N = stride = 0;
        if (dims != NULL) free(dims);
        dims = NULL;
    }

    Dataset () :dim(0), N(0), dims(NULL) {}
    Dataset (int _dim, long long _N) : dims(NULL) { reset(_dim, _N); }

    ~Dataset () { if (dims != NULL) delete[] dims; }

    /// Access the ith vector.
    const T *operator [] (int i) const {
        return (const T *)(dims + i * stride);
    }

    /// Access the ith vector.
    T *operator [] (int i) {
        return (T *)(dims + i * stride);
    }

    int getDim () const {return dim; }
    int size () const {return N; }

    void loadFvecs(const std::string &path) {


       	std::ifstream is(path.c_str(), std::ios::binary);
    	//assert(dims!=NULL);

    	char* off = dims;
    	int dimension;

    	for(long long i=0; i<this->N; i++) {
    		is.read((char*)&dimension, sizeof(dimension));
    		//assert(dimension == dim);

    		is.read(off, sizeof(T)*dim);

    		off += stride;
    	}
        std::cout << dimension << std::endl;
    	std::cout << stride << std::endl;
    	//BOOST_VERIFY(is);
    }
};

// L1 distance oracle on a dense dataset
class OracleDirect {
    const Dataset<float> &m;
public:
    OracleDirect (const Dataset<float> &m_): m(m_) {
    }
    float operator () (int i, int j) const {
        return m[i][j];
    }
};

template <typename M>
class OracleInner {
    const M &m;
};

// L1 distance oracle on a dense dataset
template <typename M>
class OracleL1 {
    const M &m;
public:
    OracleL1 (const M &m_): m(m_) {
    }
    float operator () (int i, int j) const {
        const typename M::value_type *first1 = m[i];
        const typename M::value_type *first2 = m[j];
        float r = 0.0;
        for (int i = 0; i < m.getDim(); ++i)
        {
            r += fabs(first1[i] - first2[i]);
        }
        return r;
    }
};

// L2 distance oracle on a dense dataset
// special SSE optimization is implemented for float data
template <typename M>
class OracleL2 {
    const M &m;
public:
    OracleL2 (const M &m_): m(m_) {}
    float operator () (int i, int j) const __attribute__ ((noinline));
};

template <typename M>
float OracleL2<M>::operator () (int i, int j) const {
    const typename M::value_type *first1 = m[i];
    const typename M::value_type *first2 = m[j];
    float r = 0.0;
    for (int i = 0; i < m.getDim(); ++i)
    {
        float v = first1[i] - first2[i];
        r += v * v;
    }
    return sqrt(r);
}


template <typename M>
class OracleAngular {
    const M &m;
public:
    OracleAngular (const M &m_): m(m_) {}
    float operator () (int i, int j) const __attribute__ ((noinline));
};

template <typename M>
float OracleAngular<M>::operator () (int i, int j) const {
    const typename M::value_type *first1 = m[i];
    const typename M::value_type *first2 = m[j];
    float dist_  = 0.0;
    float norm_1 = 0.0;
    float norm_2 = 0.0;
    for (int i = 0; i < m.getDim(); ++i)
    {
    	dist_  += first1[i] * first2[i];
    	norm_1 += first1[i] * first1[i];
    	norm_2 += first2[i] * first2[i];

    }
    float angle =  acos( dist_ / std::sqrt(norm_1*norm_2) );
    return angle;
}


#endif
