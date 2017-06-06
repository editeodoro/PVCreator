//--------------------------------------------------------------------
// FitsArray.hh: Definition of FitsArray class
//--------------------------------------------------------------------

#ifndef FITSARRAY_HH
#define FITSARRAY_HH

#include <iostream>
#include <string>
#include <fitsio.h>
#include <algorithm>
#include "array.hh"
#include "fitsheader.hh"

//'''''''''''''''''''''''''''''''''''''''''''''''''''''''
// FitsArray class for a generic N-Dimensional FITS file.
//_______________________________________________________

template <class T, size_t N>
class FitsArray
{	
public:

    // Constructors.
    FitsArray() {defaults();}
    FitsArray(std::string fname);
    FitsArray(size_t *dim) {setFitsArray(dim);}
    ~FitsArray() {}
    FitsArray(const FitsArray &c) {this->operator=(c);}
    void defaults() {numAxes = N; headDefined=false; verbose=true;}
    void setFitsArray(size_t *dim);

    // Overloaded operators.
    FitsArray& operator=(const FitsArray &c);
    inline T& operator() (size_t *pix) {return array(pix);}
    inline T& operator() (size_t npix) {return array(npix);}

    // Friend overloaded operators.
    template <class K, size_t M, class L> friend FitsArray<K,M>& operator+=(FitsArray<K,M>& a, const L& b);
    template <class K, size_t M, class L> friend FitsArray<K,M>& operator-=(FitsArray<K,M>& a, const L& b);
    template <class K, size_t M, class L> friend FitsArray<K,M>& operator*=(FitsArray<K,M>& a, const L& b);
    template <class K, size_t M, class L> friend FitsArray<K,M>& operator/=(FitsArray<K,M>& a, const L& b);
    template <class K, size_t M> friend FitsArray<K,M>& operator+=(FitsArray<K,M>& a, const FitsArray<K,M>& b);
    template <class K, size_t M> friend FitsArray<K,M>& operator-=(FitsArray<K,M>& a, const FitsArray<K,M>& b);
    template <class K, size_t M> friend FitsArray<K,M>& operator*=(FitsArray<K,M>& a, const FitsArray<K,M>& b);
    template <class K, size_t M> friend FitsArray<K,M>& operator/=(FitsArray<K,M>& a, const FitsArray<K,M>& b);

    /// Obvious inline functions to access a private member of class:
    int     NumAx () {return numAxes;}
    long    NumPix() {return numPix;}
    int*    AxisDim () {return axisDim;}
    int     AxesDim (int i) {return axisDim[i];}
    void    saveHead (FitsHeader &h) {head = h; headDefined=true;};
    void    setHeadDefined (bool b) {headDefined=b;}

    long    nPix  (size_t *pix) {return array.nPix(pix);}
    Array<T,N>& getArray () {return array;}
    FitsHeader& Head () {return head;}

    /// Functions for Fitsfile I/O:
    void    setArray  (const Array<T,N> input);
    void    setArray  (T *input, size_t *dim);
    bool    readArray (std::string fname);	/// Front-end to read array from Fits.
    bool    fitsread  (std::string);				  /// Read data array from Fits file.
    bool    fitswrite (std::string outfile, bool fullHead=false); /// Write a Fits Array.

	
protected:
    Array<T,N>  array;          ///< The Array data array.
    size_t      numAxes;        ///< Number of axis.
    size_t      numPix;         ///< Total number of pixel.
    size_t      axisDim[N];     ///< Array of axis dimensions of Array
    FitsHeader  head;           ///< Fits header information.
    bool        headDefined;    ///< Has been an header defined?


private:
    int     bitpix;		///< Data type code values for FITS images.
    int     datatype;	///< Data type when reading or writing data.
    bool    verbose;
};

// Overloaded operators with scalars.

template <class K, size_t M, class L>
inline FitsArray<K,M>& operator+=(FitsArray<K,M>& a, const L& b) {a.array += L(b); return a;}

template <class K, size_t M, class L>
inline FitsArray<K,M>& operator-=(FitsArray<K,M>& a, const L& b) {a.array -= L(b); return a;}

template <class K, size_t M, class L>
inline FitsArray<K,M>& operator*=(FitsArray<K,M>& a, const L& b) {a.array *= L(b); return a;}

template <class K, size_t M, class L>
inline FitsArray<K,M>& operator/=(FitsArray<K,M>& a, const L& b) {a.array /= L(b); return a;}

template <class K, size_t M, class L>
inline FitsArray<K,M> operator+(const FitsArray<K,M>& a, const L& b) {FitsArray<K,M> c=a; return (c+=b);}

template <class K, size_t M, class L>
inline FitsArray<K,M> operator-(const FitsArray<K,M>& a, const L& b) {FitsArray<K,M> c=a; return (c-=b);}

template <class K, size_t M, class L>
inline FitsArray<K,M> operator/(const FitsArray<K,M>& a, const L& b) {FitsArray<K,M> c=a; return (c*=b);}

template <class K, size_t M, class L>
inline FitsArray<K,M> operator*(const FitsArray<K,M>& a, const L& b) {FitsArray<K,M> c=a; return (c/=b);}

template <class K, size_t M, class L>
inline FitsArray<K,M> operator+(const L& b, const FitsArray<K,M>& a) {return (a+b);}

template <class K, size_t M, class L>
inline FitsArray<K,M> operator-(const L& b, const FitsArray<K,M>& a) {return (a-b);}

template <class K, size_t M, class L>
inline FitsArray<K,M> operator/(const L& b, const FitsArray<K,M>& a) {return (a*b);}

template <class K, size_t M, class L>
inline FitsArray<K,M> operator*(const L& b, const FitsArray<K,M>& a) {return (a/b);}


// Overloaded operators with other FitsArrays.

template <class K, size_t M>
inline FitsArray<K,M>& operator+=(FitsArray<K,M>& a, const FitsArray<K,M>& b) {a.array += b.array; return a;}

template <class K, size_t M>
inline FitsArray<K,M>& operator-=(FitsArray<K,M>& a, const FitsArray<K,M>& b) {a.array -= b.array; return a;}

template <class K, size_t M>
inline FitsArray<K,M>& operator*=(FitsArray<K,M>& a, const FitsArray<K,M>& b) {a.array *= b.array; return a;}

template <class K, size_t M>
inline FitsArray<K,M>& operator/=(FitsArray<K,M>& a, const FitsArray<K,M>& b) {a.array /= b.array; return a;}

template <class K, size_t M>
inline FitsArray<K,M> operator+(const FitsArray<K,M>& a, const FitsArray<K,M>& b) {FitsArray<K,M> c=a; return (c+=b);}

template <class K, size_t M>
inline FitsArray<K,M> operator-(const FitsArray<K,M>& a, const FitsArray<K,M>& b) {FitsArray<K,M> c=a; return (c-=b);}

template <class K, size_t M>
inline FitsArray<K,M> operator*(const FitsArray<K,M>& a, const FitsArray<K,M>& b) {FitsArray<K,M> c=a; return (c*=b);}

template <class K, size_t M>
inline FitsArray<K,M> operator/(const FitsArray<K,M>& a, const FitsArray<K,M>& b) {FitsArray<K,M> c=a; return (c/=b);}


//'''''''''''''''''''''''''''''''''''''''''''''''''''''''
// Derived FitsCube class for a 3-Dimensional FITS file.
//_______________________________________________________

template <class T>
class FitsCube : public FitsArray<T,3>
{
public:
    FitsCube() : FitsArray<T,3>() {}
    FitsCube(std::string fname) : FitsArray<T,3>(fname) {}
    FitsCube(size_t *d) : FitsArray<T,3>(d) {}
    FitsCube(size_t xd, size_t yd, size_t zd) {size_t d[3]={xd,yd,zd}; this->FitsCube(d);}
    FitsCube(const FitsArray<T,3> &c) {this->operator=(c);}
    ~FitsCube() {}

    inline T& operator() (size_t x, size_t y,size_t z) {return this->array(x,y,z);}
    inline T& operator() (size_t i) {return this->array(i);}
    FitsCube& operator=(const FitsArray<T,3> &c) {FitsArray<T,3>::operator =(c);}
 
    // Friend overloaded operators.
    template <class K, class L> friend FitsCube<K>& operator+=(FitsCube<K>& a, const L& b);
    template <class K, class L> friend FitsCube<K>& operator-=(FitsCube<K>& a, const L& b);
    template <class K, class L> friend FitsCube<K>& operator*=(FitsCube<K>& a, const L& b);
    template <class K, class L> friend FitsCube<K>& operator/=(FitsCube<K>& a, const L& b);
    template <class K> friend FitsCube<K>& operator+=(FitsCube<K>& a, const FitsCube<K>& b);
    template <class K> friend FitsCube<K>& operator-=(FitsCube<K>& a, const FitsCube<K>& b);
    template <class K> friend FitsCube<K>& operator*=(FitsCube<K>& a, const FitsCube<K>& b);
    template <class K> friend FitsCube<K>& operator/=(FitsCube<K>& a, const FitsCube<K>& b);

    int	DimX(){return this->axisDim[0];};
    int DimY(){return this->axisDim[1];};
    int DimZ(){return this->axisDim[2];};

    double getZphys(double z) {return (z+1-this->head.CRPIX(2))*this->head.CDELT(2)+this->head.CRVAL(2);}
    double getZgrid(double v) {return (v-this->head.CRVAL(2))/this->head.CDELT(2)+this->head.CRPIX(2)-1;}
    double* getXYZphys (double x, double y, double z) {double p[3] = {x,y,z}; double *w = new double[3]; this->head.pixToWCS(p,w); return w;}
    double* getXYphys (double x, double y) {return getXYZphys(x,y,0);}

    

};

// Overloaded operators with scalars for FitsCube class
template <class K, class L> inline FitsCube<K>& operator+=(FitsCube<K>& a, const L& b) {a.array += L(b); return a;}
template <class K, class L> inline FitsCube<K>& operator-=(FitsCube<K>& a, const L& b) {a.array -= L(b); return a;}
template <class K, class L> inline FitsCube<K>& operator*=(FitsCube<K>& a, const L& b) {a.array *= L(b); return a;}
template <class K, class L> inline FitsCube<K>& operator/=(FitsCube<K>& a, const L& b) {a.array /= L(b); return a;}
template <class K, class L> inline FitsCube<K> operator+(const FitsCube<K>& a, const L& b) {FitsCube<K> c=a; return (c+=b);}
template <class K, class L> inline FitsCube<K> operator-(const FitsCube<K>& a, const L& b) {FitsCube<K> c=a; return (c-=b);}
template <class K, class L> inline FitsCube<K> operator/(const FitsCube<K>& a, const L& b) {FitsCube<K> c=a; return (c*=b);}
template <class K, class L> inline FitsCube<K> operator*(const FitsCube<K>& a, const L& b) {FitsCube<K> c=a; return (c/=b);}
template <class K, class L> inline FitsCube<K> operator+(const L& b, const FitsCube<K>& a) {return (a+b);}
template <class K, class L> inline FitsCube<K> operator-(const L& b, const FitsCube<K>& a) {return (a-b);}
template <class K, class L> inline FitsCube<K> operator/(const L& b, const FitsCube<K>& a) {return (a*b);}
template <class K, class L> inline FitsCube<K> operator*(const L& b, const FitsCube<K>& a) {return (a/b);}
template <class K> inline FitsCube<K>& operator+=(FitsCube<K>& a, const FitsCube<K>& b) {a.array += b.array; return a;}
template <class K> inline FitsCube<K>& operator-=(FitsCube<K>& a, const FitsCube<K>& b) {a.array -= b.array; return a;}
template <class K> inline FitsCube<K>& operator*=(FitsCube<K>& a, const FitsCube<K>& b) {a.array *= b.array; return a;}
template <class K> inline FitsCube<K>& operator/=(FitsCube<K>& a, const FitsCube<K>& b) {a.array /= b.array; return a;}
template <class K> inline FitsCube<K> operator+(const FitsCube<K>& a, const FitsCube<K>& b) {FitsCube<K> c=a; return (c+=b);}
template <class K> inline FitsCube<K> operator-(const FitsCube<K>& a, const FitsCube<K>& b) {FitsCube<K> c=a; return (c-=b);}
template <class K> inline FitsCube<K> operator*(const FitsCube<K>& a, const FitsCube<K>& b) {FitsCube<K> c=a; return (c*=b);}
template <class K> inline FitsCube<K> operator/(const FitsCube<K>& a, const FitsCube<K>& b) {FitsCube<K> c=a; return (c/=b);}


//'''''''''''''''''''''''''''''''''''''''''''''''''''''''
// Derived FitsImage class for a 2-Dimensional FITS file.
//_______________________________________________________

template <class T>
class FitsImage : public FitsArray<T,2>
{
public:
    FitsImage() : FitsArray<T,2>() {}
    FitsImage(std::string fname) : FitsArray<T,2>(fname) {}
    FitsImage(size_t xd, size_t yd) {size_t d[2]={xd,yd}; FitsArray<T,2>::FitsArray(d);}
    FitsImage(size_t *d) : FitsArray<T,2>::FitsArray(d) {}
    FitsImage(const FitsArray<T,2> &c) {this->operator=(c);}
    ~FitsImage() {}

    inline T& operator() (size_t x, size_t y) {return this->array(x,y);};
    FitsImage& operator=(const FitsArray<T,2> &c) {FitsArray<T,2>::operator=(c);}

    // Friend overloaded operators.
    template <class K, class L> friend FitsImage<K>& operator+=(FitsImage<K>& a, const L& b);
    template <class K, class L> friend FitsImage<K>& operator-=(FitsImage<K>& a, const L& b);
    template <class K, class L> friend FitsImage<K>& operator*=(FitsImage<K>& a, const L& b);
    template <class K, class L> friend FitsImage<K>& operator/=(FitsImage<K>& a, const L& b);
    template <class K> friend FitsImage<K>& operator+=(FitsImage<K>& a, const FitsImage<K>& b);
    template <class K> friend FitsImage<K>& operator-=(FitsImage<K>& a, const FitsImage<K>& b);
    template <class K> friend FitsImage<K>& operator*=(FitsImage<K>& a, const FitsImage<K>& b);
    template <class K> friend FitsImage<K>& operator/=(FitsImage<K>& a, const FitsImage<K>& b);
    
    
    int DimX(){return this->axisDim[0];}
    int DimY(){return this->axisDim[1];}

    double getXphys(double x) {return (x+1-this->head.CRPIX(0))*this->head.CDELT(0)+this->head.CRVAL(0);}
    double getYphys(double y) {return (y+1-this->head.CRPIX(1))*this->head.CDELT(1)+this->head.CRVAL(1);}

};

// Overloaded operators with scalars for FitsImage class
template <class K, class L> inline FitsImage<K>& operator+=(FitsImage<K>& a, const L& b) {a.array += L(b); return a;}
template <class K, class L> inline FitsImage<K>& operator-=(FitsImage<K>& a, const L& b) {a.array -= L(b); return a;}
template <class K, class L> inline FitsImage<K>& operator*=(FitsImage<K>& a, const L& b) {a.array *= L(b); return a;}
template <class K, class L> inline FitsImage<K>& operator/=(FitsImage<K>& a, const L& b) {a.array /= L(b); return a;}
template <class K, class L> inline FitsImage<K> operator+(const FitsImage<K>& a, const L& b) {FitsImage<K> c=a; return (c+=b);}
template <class K, class L> inline FitsImage<K> operator-(const FitsImage<K>& a, const L& b) {FitsImage<K> c=a; return (c-=b);}
template <class K, class L> inline FitsImage<K> operator/(const FitsImage<K>& a, const L& b) {FitsImage<K> c=a; return (c*=b);}
template <class K, class L> inline FitsImage<K> operator*(const FitsImage<K>& a, const L& b) {FitsImage<K> c=a; return (c/=b);}
template <class K, class L> inline FitsImage<K> operator+(const L& b, const FitsImage<K>& a) {return (a+b);}
template <class K, class L> inline FitsImage<K> operator-(const L& b, const FitsImage<K>& a) {return (a-b);}
template <class K, class L> inline FitsImage<K> operator/(const L& b, const FitsImage<K>& a) {return (a*b);}
template <class K, class L> inline FitsImage<K> operator*(const L& b, const FitsImage<K>& a) {return (a/b);}
template <class K> inline FitsImage<K>& operator+=(FitsImage<K>& a, const FitsImage<K>& b) {a.array += b.array; return a;}
template <class K> inline FitsImage<K>& operator-=(FitsImage<K>& a, const FitsImage<K>& b) {a.array -= b.array; return a;}
template <class K> inline FitsImage<K>& operator*=(FitsImage<K>& a, const FitsImage<K>& b) {a.array *= b.array; return a;}
template <class K> inline FitsImage<K>& operator/=(FitsImage<K>& a, const FitsImage<K>& b) {a.array /= b.array; return a;}
template <class K> inline FitsImage<K> operator+(const FitsImage<K>& a, const FitsImage<K>& b) {FitsImage<K> c=a; return (c+=b);}
template <class K> inline FitsImage<K> operator-(const FitsImage<K>& a, const FitsImage<K>& b) {FitsImage<K> c=a; return (c-=b);}
template <class K> inline FitsImage<K> operator*(const FitsImage<K>& a, const FitsImage<K>& b) {FitsImage<K> c=a; return (c*=b);}
template <class K> inline FitsImage<K> operator/(const FitsImage<K>& a, const FitsImage<K>& b) {FitsImage<K> c=a; return (c/=b);}



template <class T> inline int selectBITPIX();
template <> inline int selectBITPIX<short>() {return SHORT_IMG;};
template <> inline int selectBITPIX<int>() {return SHORT_IMG;};
template <> inline int selectBITPIX<long>() {return LONG_IMG;};
template <> inline int selectBITPIX<float>() {return FLOAT_IMG;};
template <> inline int selectBITPIX<double>() {return DOUBLE_IMG;};


template <class T> inline int selectDATATYPE();
template <> inline int selectDATATYPE<short>() {return TSHORT;};
template <> inline int selectDATATYPE<int>() {return TINT;};
template <> inline int selectDATATYPE<long>() {return TLONG;};
template <> inline int selectDATATYPE<float>() {return TFLOAT;};
template <> inline int selectDATATYPE<double>() {return TDOUBLE;};


template <class T, size_t N>
inline FitsArray<T,N>::FitsArray(std::string fname) {
    defaults();
    if (!readArray(fname)) std::cerr << "\nFitsArray: cannot read " << fname << std::endl;
}


template <class T, size_t N>
inline FitsArray<T,N>& FitsArray<T,N>::operator=(const FitsArray<T,N> &c) {
    
    if(this==&c) return *this;
    numPix	= c.numPix;
    numAxes	= c.numAxes;
    for (int i=0; i<numAxes;  i++) axisDim[i] = c.axisDim[i];
    array=c.array;
    headDefined = c.headDefined;
    if (headDefined) head = c.head;
    return *this;
}


template <class T, size_t N>
inline void FitsArray<T,N>::setFitsArray(size_t *dim) {
    defaults();
    numPix = 1;
    for (int i=N; i--;) {
        axisDim[i]=dim[i];
        numPix *= axisDim[i];
    }
    array.setArray(dim);
    head.setNumAx(N);
}

/**=====================================================================================*/ 
/** FUNCTIONS FOR INPUT AND OUTPUT */

template <class T, size_t N>
inline void FitsArray<T,N>::setArray (const Array<T,N> input) {

    numAxes = N;
    for (int i=N; i--;) axisDim[i]=input.Dim(i);
    numPix = input.Size();
    array = input;
}



template <class T, size_t N>
void FitsArray<T,N>::setArray (T *input, size_t *dim) {
    Array<T,N> d(dim, input);
    setArray(d);
}



template <class T, size_t N>
bool FitsArray<T,N>::readArray (std::string fname) {

    if(!head.header_read(fname)) return false;
    headDefined = true;
    if (head.NumAx()>N) {
        std::cout << "\nFitsArray WARNING: Reading " << N
                  << " out of " << head.NumAx() << " dimensions. \n";
    }
    else if (head.NumAx()<N) {
        std::cout << "\nFitsArray ERROR: Cannot read " << head.NumAx()
                  << "D object into a " << N << "D object. \n";
        return false;
    }
    numPix=1;
    for (int i=0; i<numAxes; i++) {
        axisDim[i] = head.NAXIS(i);
        numPix *= axisDim[i];
    }
    if (!fitsread(fname)) return false;
    return true;
}



template <class T, size_t N>
bool FitsArray<T,N>::fitsread(std::string fname) {

    fitsfile *fptr;
    int status=0, anynul, fpixel=1;
    float nulval;

    if (verbose) {
        std::cout << "\nOpening file "<< fname << std::endl;
        std::cout << "Reading ";
        for (int i=0; i<numAxes; i++) {
            std::cout << axisDim[i];
            if (i!=numAxes-1) std::cout << " x ";
        }
        std::cout << " pixels FITS file... ";
    }

    // Open the FITS file
    if(fits_open_file(&fptr, fname.c_str(), READONLY, &status) ){
        fits_report_error(stderr, status);
        return false;
    }

    // Read elements from the FITS data array
    if (!array.P()) array.setArray(axisDim);

    status=0;
    if (fits_read_img(fptr, selectDATATYPE<T>(), fpixel, numPix, &nulval, array.P(), &anynul, &status)){
        fits_report_error(stderr, status);
        return false;
    }

    // Close the FITS File
    if (fits_close_file(fptr, &status)) fits_report_error(stderr, status);


    if (verbose) std::cout << "Done.\n" << std::endl;

    return true;
}



template <class T, size_t N>
bool FitsArray<T,N>::fitswrite(std::string outfile, bool fullHead) {

    fitsfile *fptr;
    int fpixel = 1, status=0;    
    long dnaxes[numAxes];
    for (int i=numAxes; i--;) dnaxes[i]=axisDim[i];

    remove(outfile.c_str());

    if (fits_create_file(&fptr, outfile.c_str(), &status)) {
        fits_report_error(stderr, status);
        return false;
    }

    status=0;
    if (fits_create_img(fptr, selectBITPIX<T>(), numAxes, dnaxes, &status)) {
        fits_report_error(stderr, status);
        return false;
    }
    
    if (headDefined) head.header_write(fptr,fullHead);

    status=0;
    if (fits_write_img(fptr, selectDATATYPE<T>(), fpixel, numPix, array.P(), &status)) {
        fits_report_error(stderr, status);
        return false;
    }

    if (fits_close_file(fptr, &status)) {
        fits_report_error(stderr, status);
    }

    return true;
}

#endif
