//-----------------------------------------------------------------------
// fitsheader.hh: Definition of FitsHeader class, a class to collect Header info.
//-----------------------------------------------------------------------

#ifndef FITSHEADER_HH
#define FITSHEADER_HH

#include <iostream>
#include <cstring>
#include <cmath>
#include <vector>
#include <fitsio.h>
#include <cctype>
#include <iomanip>
#include <sstream>
#include <type_traits>
#include <wcslib/wcs.h>
#include <wcslib/wcsunits.h>
#include <wcslib/wcshdr.h>
#include <wcslib/wcsfix.h>

//enum CUNITS {DEGREE,ASEC,AMIN,M_S,KM_S,HZ,KHZ,MHZ,GHZ,MUM};
//enum CTYPES {RA,DEC,LAT,LONG,VELO,FREQ,WAVE};

class FitsHeader
{
public: 
    FitsHeader();                               /// Default constructor.
    FitsHeader(int numaxis) {setNumAx(numaxis);}
    virtual ~FitsHeader();                       /// Default destructor.
    FitsHeader(const FitsHeader& h);		/// Copy constructor.
    FitsHeader& operator= (const FitsHeader& h);/// Assignement operator.

    double operator()(std::string s);
    double operator[](std::string s) {this->operator()(s);}

    /// Obvious inline functions
    int     NumAx () {return numAxes;};
    int     BITPIX () {return bitpix;};
    long    *DIMAXES () {return dimAxes;};
    long&   NAXIS (int i) {return dimAxes[i];};
    double  *CRPIX () {return crpix;};
    double  *CRVAL () {return crval;};
    double  *CDELT () {return cdelt;};
    double  BMAJ () {return bmaj;};
    double  BMIN () {return bmin;};
    double  BPA	 () {return bpa;};
    float   BeamArea() {return beamArea;};
    float   BZERO () {return bzero;};
    float   BSCALE () {return bscale;};
    float   BLANK () {return blank;};
    float   EPOCH () {return epoch;};
    double  FREQ0 () {return freq0;};
    double  CROTA () {return crota;};
    double  DATAMAX () {return datamax;};
    double  DATAMIN () {return datamin;};
    double  CDELT (int i) {return cdelt[i];};
    double  CRVAL (int i) {return crval[i];};
    double  CRPIX (int i) {return crpix[i];};
    double  DRVAL3 () {return drval3;};
    double  PixScale () {return (fabs(cdelt[0])+fabs(cdelt[1]))/2.;};
    struct  wcsprm *WCS () {return wcs;};

    std::vector<std::string>& Keys () {std::vector<std::string> &k=keys; return k;};
    std::string OBJECT () {return object;};
    std::string BUNIT () {return bunit;};
    std::string BTYPE () {return btype;};
    std::string CTYPE (int i) {return ctype[i];};
    std::string CUNIT (int i) {return cunit[i];};
    std::string DUNIT3 () {return dunit3;};
    std::string TELESCOPE () {return telescope;};

    void    setBitpix (int i) {bitpix = i;};
    void    setDimAx (int i, long val) {dimAxes[i] = val;};
    void    setCrpix (int i, float val) {crpix[i] = val;};
    void    setCrval (int i, float val) {crval[i] = val;};
    void    setCdelt (int i, float val) {cdelt[i] = val;};
    void    setDrval3 (double val) {drval3=val;};
    void    setDunit3 (std::string s) {dunit3=s;};
    void    setBmaj  (float val) {bmaj = val;};
    void    setBmin  (float val) {bmin = val;};
    void    setBpa   (float val) {bpa = val;};
    void    setBzero (float val) {bzero = val;};
    void    setBscale(float val) {bscale = val;};
    void    setBlank (float val) {blank = val;};
    void    setEpoch (float val) {epoch = val;};
    void    setBunit (std::string ch) {bunit = ch;};
    void    setDataMax (double val) {datamax=val;};
    void    setDataMin (double val) {datamin=val;};
    void    setMinMax (double minn, double maxx) {datamin=minn;datamax=maxx;};
    void    setFreq0 (double val) {freq0=val;};
    void    setCtype (int i, std::string s) {ctype[i] = s;};
    void    setCunit (int i, std::string s) {cunit[i] = s;};
    void    setBtype (std::string s) {btype = s;};
    void    setName  (std::string s) {object = s;};
    void    setTelesc(std::string s) {telescope = s;};
    void    setWarning (bool b) {warning=b;};

    void    Warning(std::string s) {if (warning) std::cout << s << std::endl;};


    /// Functions defined in FitsHeader.cpp.
    void    setNumAx (int n);
    void    calcArea ();                                       /// Calculate beam area from bmaj & bmin.
    bool    header_read (std::string fname);                   /// Read from FitsHeader of a Fits file.
    void    header_write (fitsfile *fptr, bool fullHead);      /// Write FitsHeader of a Fits FIle.
    int     wcsToPix(const double *world, double *pix, size_t npts=1);  /// Convert WCS to pixel
    int     pixToWCS(const double *pix, double *world, size_t npts=1);  /// Convert a n-dim pixel to WCS
        
    template <class T>							/// Read the request keyword and write on "key".
    bool read_keyword(std::string keyword, T &key, bool err=false);

private:
    int     numAxes;        ///< Number of axes.
    int     bitpix;         ///< Image type.
    long    *dimAxes=NULL;  ///< Dimensions of axes.
    double  *crpix=NULL;    ///< Central pixels.
    double  *crval=NULL;    ///< Values of central pixels.
    double  *cdelt=NULL;    ///< Delta pixel.
    double  drval3;         ///< Secondary reference value of third axis.
    double  bmaj;           ///< The major main-beam FWHM.
    double  bmin;           ///< The minor main-beam FWHM.
    double  bpa;            ///< The beam position angle.
    float   bzero;          ///< Bias for real values.
    float   bscale;         ///< Scale for physical values.
    float   blank;          ///< Value for blank pixel.
    float   beamArea;       ///< The area of the beam.
    float   epoch;          ///< Epoch for coordinates.
    double  freq0;          ///< Frequency at rest.
    double  crota;          ///< Rotation angle.
    double  datamin;        ///< Minimum pixel value.
    double      datamax;    ///< Maximum data value.
    std::string fitsname;   ///< The name of the fitsfile.
    std::string	btype;      ///< Beam type.
    std::string bunit;      ///< Units of pixel value.
    std::string object;     ///< The name of the object.
    std::string *ctype;     ///< Type of axis.
    std::string *cunit;     ///< Unity of axis.
    std::string dunit3;     ///< Secondary units of third axis.
    std::string telescope;  ///< Instrument.
    std::vector<std::string> keys;	///< Whole FitsHeader as strings.

    struct wcsprm *wcs;     ///< The WCS parameters in a struct from the wcslib library.
    int    nwcs;            ///< The number of WCS parameters
    bool   wcsIsGood;       ///< A flag indicating whether there is a valid WCS
    
    bool    warning;               ///< Write warning on std::cout.
};


inline FitsHeader::FitsHeader () {

    bitpix = FLOAT_IMG;
    numAxes = bmaj = bmin = bpa = beamArea = freq0 = 0.;
    datamin = datamax = 0.;
    dunit3 = "";
    object = "NONE";
    warning = true;
    wcs = new struct wcsprm;
    wcs->flag=-1;
    wcsIsGood=false;
    nwcs=0;
}


inline FitsHeader::~FitsHeader () {

    if (dimAxes) {
        delete [] dimAxes;
        delete [] crpix;
        delete [] crval;
        delete [] cdelt;
        delete [] ctype;
        delete [] cunit;
    }

    keys.clear();
    wcsvfree(&nwcs,&wcs);
}


inline FitsHeader::FitsHeader(const FitsHeader& h) {

    operator=(h);
}


inline FitsHeader& FitsHeader::operator=(const FitsHeader& h) {

    if(this == &h) return *this;
    this->numAxes = h.numAxes;
    this->bitpix  = h.bitpix;

    if (dimAxes) {
        delete [] this->dimAxes;
        delete [] this->crpix;
        delete [] this->crval;
        delete [] this->cdelt;
        delete [] this->ctype;
        delete [] this->cunit;
    }

    if(h.dimAxes) {
        this->dimAxes = new long[numAxes];
        this->crpix	= new double[numAxes];
        this->crval	= new double[numAxes];
        this->cdelt	= new double[numAxes];
        this->cunit = new std::string[numAxes];
        this->ctype = new std::string[numAxes];
        for (int i=0; i<numAxes; i++) {
            this->dimAxes[i] = h.dimAxes[i];
            this->crpix[i] = h.crpix[i];
            this->crval[i] = h.crval[i];
            this->cdelt[i] = h.cdelt[i];
            this->ctype[i] = h.ctype[i];
            this->cunit[i] = h.cunit[i];
        }
    }

    this->bmaj      = h.bmaj;
    this->bmin      = h.bmin;
    this->bpa       = h.bpa;
    this->bzero     = h.bzero;
    this->bscale    = h.bscale;
    this->blank     = h.blank;
    this->beamArea  = h.beamArea;
    this->epoch     = h.epoch;
    this->freq0     = h.freq0;
    this->fitsname  = h.fitsname;
    this->bunit     = h.bunit;
    this->btype     = h.btype;
    this->object    = h.object;
    this->telescope = h.telescope;
    this->dunit3    = h.dunit3;
    this->drval3    = h.drval3;
    this->datamin   = h.datamin;
    this->datamax   = h.datamax;
    this->warning   = h.warning;


    this->wcs = new struct wcsprm;
    this->wcs->flag = -1;
    wcsini(true, h.wcs->naxis, this->wcs);
    wcscopy(true, h.wcs, this->wcs);
    wcsset(this->wcs);
    this->nwcs      = h.nwcs;
    this->wcsIsGood = h.wcsIsGood;

    for (unsigned int i=0; i<h.keys.size(); i++)
        this->keys.push_back(h.keys[i]);

    calcArea();

    return *this;

}


inline double FitsHeader::operator()(std::string s) {

    if (s=="BITPIX") return bitpix;
    if (s=="NAXIS" ) return numAxes;
    if (s=="NAXIS1" && numAxes>0) return dimAxes[0];
    if (s=="NAXIS2" && numAxes>1) return dimAxes[1];
    if (s=="NAXIS3" && numAxes>2) return dimAxes[2];
    if (s=="NAXIS4" && numAxes>3) return dimAxes[3];
    if (s=="CRPIX1" && numAxes>0) return crpix[0];
    if (s=="CRPIX2" && numAxes>1) return crpix[1];
    if (s=="CRPIX3" && numAxes>2) return crpix[2];
    if (s=="CRPIX4" && numAxes>3) return crpix[3];
    if (s=="CRVAL1" && numAxes>0) return crval[0];
    if (s=="CRVAL2" && numAxes>1) return crval[1];
    if (s=="CRVAL3" && numAxes>2) return crval[2];
    if (s=="CRVAL4" && numAxes>3) return crval[3];
    if (s=="CDELT1" && numAxes>0) return cdelt[0];
    if (s=="CDELT2" && numAxes>1) return cdelt[1];
    if (s=="CDELT3" && numAxes>2) return cdelt[2];
    if (s=="CDELT4" && numAxes>3) return cdelt[3];
    if (s=="BMAJ"||s=="BMMAJ"||s=="CLBMMAJ") return bmaj;
    if (s=="BMIN"||s=="BMMIN"||s=="CLBMMIN") return bmin;
    if (s=="BPA"||s=="BMPA"||s=="CLBMPA") return bpa;
    if (s=="FREQ0"||s=="RESTFREQ"||s=="RESTFRQ") return freq0;
    if (s=="BZERO")  return bzero;
    if (s=="BSCALE") return bscale;
    if (s=="BLANK")  return blank;
    if (s=="EPOCH")  return epoch;
    if (s=="CROTA")  return crota;
    if (s=="DATAMAX") return datamax;
    if (s=="DATAMIN") return datamin;

    double val=0;
    if(read_keyword(s,val)) return val;
    else {
        std::cerr << "\n FitsHeader ERROR: cannot read key " << s << std::endl;
        return 0;
    }
}


inline void FitsHeader::setNumAx (int size){

    if (dimAxes) {
        delete [] dimAxes;
        delete [] crpix;
        delete [] crval;
        delete [] cdelt;
        delete [] cunit;
        delete [] ctype;
    }
    numAxes = size;
    dimAxes = new long[numAxes];
    crpix   = new double[numAxes];
    crval   = new double[numAxes];
    cdelt   = new double[numAxes];
    cunit   = new std::string[numAxes];
    ctype   = new std::string[numAxes];

}


inline void FitsHeader::calcArea () {

    float AvPixScale = (fabs(cdelt[0])+fabs(cdelt[1]))/2.;
    float scalbmaj = bmaj/AvPixScale;
    float scalbmin = bmin/AvPixScale;
    beamArea =  M_PI * (scalbmaj/2.) * (scalbmin/2.) / M_LN2;

}


inline bool FitsHeader::header_read (std::string fname) {

    fitsfile *fptr;
    int status=0, nfound;
    char comment[72];


    fitsname = fname;
    char filename[200];
    strcpy(filename, fname.c_str());

    if (fits_open_file(&fptr, filename, READONLY, &status)) {
        fits_report_error(stderr, status);
        return false;
    }

    status=0;
    if (fits_get_img_type(fptr, &bitpix, &status)) {
        fits_report_error(stderr, status);
        bitpix = FLOAT_IMG;
    }

    status=0;
    if (fits_get_img_dim(fptr, &numAxes, &status))
        fits_report_error(stderr, status);

    if (!dimAxes) {
        dimAxes = new long [numAxes];
        crpix = new double [numAxes];
        crval = new double [numAxes];
        cdelt = new double [numAxes];
        ctype = new std::string[numAxes];
        cunit = new std::string[numAxes];
    }

    status = 0;
    if(fits_get_img_size(fptr, numAxes, dimAxes, &status)){
        fits_report_error(stderr, status);
    }

    char Bunit[20], Btype[20], name[20], Tel[20], Dunit3[20], Keys[100];
    for (int i=0; i<20; i++) {
        Bunit[i]=Btype[i]=name[i]=Tel[i]=Dunit3[i]=Keys[i]=' ';
    }
    for (int i=20; i<100; i++) Keys[i]=' ';

    int nkeys;
    fits_get_hdrspace(fptr, &nkeys, NULL, &status);
    for (int i=1; i<=nkeys; i++) {
        fits_read_record(fptr,i,Keys,&status);
        keys.push_back(Keys);
    }

    status=0;
    fits_read_keys_dbl (fptr, "CDELT", 1, numAxes, cdelt, &nfound, &status);
    if (nfound==0) fits_report_error(stderr, status);

    status=0;
    fits_read_keys_dbl (fptr, "CRPIX", 1, numAxes, crpix, &nfound, &status);
    if (nfound==0) fits_report_error(stderr, status);

    status=0;
    fits_read_keys_dbl (fptr, "CRVAL", 1, numAxes, crval, &nfound, &status);
    if (nfound==0) fits_report_error(stderr, status);

    status=0;
    if (fits_read_key_dbl (fptr, "DRVAL3", &drval3, comment, &status)) {
        if (status!=202) fits_report_error(stderr, status);
        drval3=0;
    }

    char **Ctype = new char*[numAxes];
    char **Cunit = new char*[numAxes];
    for (int i=0; i<numAxes; i++) {
        Ctype[i] = new char[25];
        Cunit[i] = new char[25];
    }

    status=0;
    fits_read_keys_str (fptr, "CTYPE", 1, numAxes, Ctype, &nfound, &status);
    if (nfound==0) {
        Warning("Error reading FitsHeader (CTYPEs). Assuming [RA,DEC,VELO].");
        if (numAxes>0) ctype[0] = "RA---NCP";
        if (numAxes>1) ctype[1] = "DEC--NCP";
        if (numAxes>2) ctype[2] = "VELO-HELO";
    }
    else {
        for (int i=0; i<numAxes; i++) {
            ctype[i] = Ctype[i];
        }
    }

    if (ctype[0].find("RA")>=0 && crval[0]<0) crval[0]+=360.;

    status=0;
    fits_read_keys_str (fptr, "CUNIT", 1, numAxes, Cunit, &nfound, &status);
    if (nfound==0) {
        if (numAxes>0) cunit[0] = "DEGREE";
        if (numAxes>1) cunit[1] = "DEGREE";
        if (numAxes>2) {
            if (ctype[2]=="FREQ" || ctype[2]=="freq" || ctype[2]=="Freq") {
                cunit[2] = "HZ";
                Warning("Error reading FitsHeader (CUNITs). Assuming [DEGREE,DEGREE,HZ]");
            }
            else {
                cunit[2] = "M/S";
                Warning("Error reading FitsHeader (CUNITs). Assuming [DEGREE,DEGREE,M/S]");
            }
        }
    }
    else {
        for (int i=0; i<numAxes; i++) {
            cunit[i] = Cunit[i];
        }
    }

    for (int i=0; i<numAxes; i++) {
        delete [] Ctype[i];
        delete [] Cunit[i];
    }
    delete [] Cunit;
    delete [] Ctype;

    status=0;
    if (fits_read_key_str (fptr, "DUNIT3", Dunit3, comment, &status)) {
        if (status!=202) fits_report_error(stderr, status);
        dunit3 = "NONE";
    }
    else dunit3 = Dunit3;

    status=0;
    if (fits_read_key_str (fptr, "BUNIT", Bunit, comment, &status)) {
        if (status==202)
            Warning("Error reading FitsHeader (BUNIT): keyword not found.");
        else fits_report_error(stderr, status);
        bunit = "NONE";
    }
    else bunit = Bunit;

    status=0;
    if (fits_read_key_str (fptr, "BTYPE", Btype, comment, &status)) {
        btype = "NONE";
    }
    else btype = Btype;

    status=0;
    if (fits_read_key_flt (fptr, "BZERO", &bzero, comment, &status)) {
        bzero = 0;
    }

    status=0;
    if (fits_read_key_flt (fptr, "BSCALE", &bscale, comment, &status)) {
        bscale = 1;
    }

    status=0;
    if (fits_read_key_flt (fptr, "BLANK", &blank, comment, &status)) {
        blank = 0;
    }

    status=0;
    if (fits_read_key_flt (fptr, "EPOCH", &epoch, comment, &status)) {
        if (fits_read_key_flt (fptr, "EQUINOX", &epoch, comment, &status))
            epoch = 0;
    }

    status=0;
    if (fits_read_key_dbl (fptr, "DATAMIN", &datamin, comment, &status)) {
        datamin = 0;
    }

    status=0;
    if (fits_read_key_dbl (fptr, "DATAMAX", &datamax, comment, &status)) {
        datamax = 0;
    }

    status=0;
    if (fits_read_key_dbl (fptr, "FREQ0", &freq0, comment, &status)) {
        status=0;
        if (fits_read_key_dbl (fptr, "RESTFREQ", &freq0, comment, &status)) {
            status=0;
            if (fits_read_key_dbl (fptr, "RESTFRQ", &freq0, comment, &status)) {
                if (dunit3=="NONE" || drval3==0 || (cunit[2]!="HZ" && cunit[2]!="hz" &&
                                                    cunit[2]!="MHZ" && cunit[2]!="Mhz" && cunit[2]!="mhz" &&
                                                    cunit[2]!="GHZ" && cunit[2]!="Ghz" && cunit[2]!="ghz" &&
                                                    dunit3!="KM/S" && dunit3!="km/s" && dunit3!="M/S"  && dunit3!="m/s")) {
                    Warning("Error reading FitsHeader (FREQ0-RESTFREQ) Assuming 1.4204057 GHz.");
                    freq0 = 0.1420405751786E10;
                }
                else {
                    double drval3ms=0.;
                    double crval3hz=0.;
                    if (dunit3=="KM/S" || dunit3=="km/s") drval3ms=drval3*1000;
                    else if (dunit3=="M/S" || dunit3=="m/s") drval3ms=drval3;
                    if (cunit[2]=="HZ" || cunit[2]=="hz") crval3hz = crval[2];
                    else if (cunit[2]=="MHZ" || cunit[2]=="Mhz" || cunit[2]=="mhz") crval3hz=crval[2]*1.E06;
                    else if (cunit[2]=="GHZ" || cunit[2]=="Ghz" || cunit[2]=="ghz") crval3hz=crval[2]*1.E09;
                    freq0 = crval3hz*sqrt((299792458.+drval3ms)/(299792458.-drval3ms));
                }
            }
        }
    }

    status=0;
    if (fits_read_key_dbl (fptr, "CROTA", &crota, comment, &status)) {
        crota = 0;
    }

    status=0;
    if (fits_read_key_str (fptr, "OBJECT", name, comment, &status)) {
        if (status==202) Warning("Error reading FitsHeader (OBJECT): keyword not found.");
        else fits_report_error(stderr, status);
        object = "NONE";
    }
    else object = name;

    for (uint i=0; i<object.size();i++) {
        if (object[i]=='/' || object[i]=='\\') object.replace(i,1, "-");
        if (isspace(object[i])) object.replace(i,1, "_");
    }

    status = 0;
    if (fits_read_key_str (fptr, "TELESCOP", Tel, comment, &status)) {
        status = 0;
        if (fits_read_key_str (fptr, "INSTRUME", Tel, comment, &status)) {
            if (status==202) Warning("Error reading FitsHeader (TELESCOP-INSTRUME): keyword not found.");
            else fits_report_error(stderr, status);
            telescope = "NONE";
        }
        else telescope = Tel;
    }
    else telescope = Tel;


    double clbmaj=0, clbmin=0;
    status=0;
    if (fits_read_key_dbl (fptr, "BMAJ", &bmaj, comment, &status)) {
        status = 0;
        if (fits_read_key_dbl (fptr, "BMMAJ", &bmaj, comment, &status)) {
            status = 0;
            if (fits_read_key_dbl (fptr, "CLBMMAJ", &clbmaj, comment, &status)) {
                if (status==202) Warning("Error reading FitsHeader (BMAJ-BMMAJ-CLBMMAJ): keyword not found.");
                else fits_report_error(stderr, status);
                bmaj = 0;
            }
        }
    }

    status=0;
    if (fits_read_key_dbl (fptr, "BMIN", &bmin, comment, &status)) {
        status = 0;
        if (fits_read_key_dbl (fptr, "BMMIN", &bmin, comment, &status)) {
            status = 0;
            if (fits_read_key_dbl (fptr, "CLBMMIN", &clbmin, comment, &status)) {
                if (status==202) Warning("Error reading FitsHeader (BMIN-BMMIN-CLBMMIN): keyword not found.");
                else fits_report_error(stderr, status);
                bmin = 0;
            }
        }
    }

    status=0;
    if (fits_read_key_dbl (fptr, "BPA", &bpa, comment, &status)) {
        status = 0;
        if (fits_read_key_dbl (fptr, "BMPA", &bpa, comment, &status)) {
            if (status==202) Warning("Error reading FitsHeader (BPA-BMPA): keyword not found.");
            else fits_report_error(stderr, status);
            bpa = 0;
        }
    }

    if (bmaj==0 && bmin==0 && bpa==0 && clbmaj==0 && clbmin==0) {
        for (unsigned int i=0; i<keys.size(); i++) {
            int found = keys[i].find("BMAJ=");
            char *pEnd;
            if (found>=0) {
                pEnd = &keys[i].at(found+5);
                bmaj = strtod(pEnd,NULL);
            }
            found = keys[i].find("BMIN=");
            if (found>=0) {
                pEnd = &keys[i].at(found+5);
                bmin = strtod(pEnd,NULL);

            }
            found = keys[i].find("BPA=");
            if (found>=0) {
                pEnd = &keys[i].at(found+4);
                bpa = strtod(pEnd,NULL);
            }
        }
        if (bmaj!=0 && bmin!=0) {
            std::cout << "\n--------> WARNING: beam information found in HISTORY keywords: <--------\n"
                      << " BMAJ = " << bmaj << " " << cunit[0]
                      << "  BMIN = " << bmin << " " << cunit[0]
                      << "  BPA = "  << bpa  << " DEGREE\n";
            std::cout << " It is heartly recommended to check these values before going on!!\n\n";
        }
    }

    char bm[30];
    status=0;
    fits_read_key_str (fptr, "BMMAJ", bm, comment, &status);
    std::string bmstr = bm;
    bool arcsecbeam = false;
    int found = bmstr.find("D");
    if (found>=0) arcsecbeam = true;
    if (arcsecbeam) bmaj /= 3600.;

    status=0;
    fits_read_key_str (fptr, "BMMIN", bm, comment, &status);
    bmstr = bm;
    arcsecbeam = false;
    found = bmstr.find("D");
    if (found>=0) arcsecbeam = true;
    if (arcsecbeam) bmin /= 3600.;

    if (clbmaj!=0) bmaj = clbmaj/3600.;
    if (clbmin!=0) bmin = clbmin/3600.;



    int noComments = 1;     // fits_hdr2str will ignore COMMENT, HISTORY etc
    int nExc = 0;
    char *hdr=0;

    // Read in the entire PHU of the FITS file to a std::string.
    // This will be read by the wcslib functions to extract the WCS.
    status = 0;
    fits_hdr2str(fptr, noComments, NULL, nExc, &hdr, &nkeys, &status);

    status = wcsini(true, numAxes, wcs);

    int relax=1; // for wcspih -- admit all recognised informal WCS extensions
    int ctrl=2;  // for wcspih -- report each rejected card and its reason for rejection
    int nreject;
    // Parse the FITS FitsHeader to fill in the wcsprm structure
    status=wcspih(hdr, nkeys, relax, ctrl, &nreject, &nwcs, &wcs);


    int stat[NWCSFIX];
    // Applies all necessary corrections to the wcsprm structure
    //  (missing cards, non-standard units or spectral types, ...)
    status = wcsfix(1, (const int*)dimAxes, wcs, stat);

    // Set up the wcsprm struct. Report if something goes wrong.
    status = wcsset(wcs);
    // Re-do the corrections to account for things like NCP projections
    status = wcsfix(1, (const int*)dimAxes, wcs, stat);

    char stype[5],scode[5],sname[22],units[8],ptype,xtype;
    int restreq;

    status = spctyp(wcs->ctype[wcs->spec],stype,scode,sname,units,&ptype,&xtype,&restreq);


    // Close the FITS File
    status=0;
    if (fits_close_file(fptr, &status))
        fits_report_error(stderr, status);

    calcArea();

    return true;

}


inline void FitsHeader::header_write (fitsfile *fptr, bool fullHead) {

    int status=0;
    char com[]= "  ";

    for (int i=0; i<numAxes; i++) {
        std::stringstream ss;
        ss << std::fixed << std::setprecision(0) << i+1;
        std::string key_towrite[5] = {"CRPIX","CRVAL","CDELT","CTYPE","CUNIT"};
        double pp[3] = {crpix[i],crval[i],cdelt[i]};
        std::string p[2] = {ctype[i],cunit[i]};
        for (int j=0; j<3; j++) {
            key_towrite[j] += ss.str();
            fits_update_key_dbl(fptr, key_towrite[j].c_str(), pp[j], 10, com, &status);
        }
        for (int j=3; j<5; j++) {
            key_towrite[j] += ss.str();
            fits_update_key_str(fptr, key_towrite[j].c_str(), p[j-3].c_str(), com, &status);
        }
    }

    if (drval3!=0) fits_update_key_dbl(fptr, "DRVAL3", drval3, 10, com, &status);
    if (dunit3!="NONE" && dunit3!="") fits_update_key_str(fptr, "DUNIT3", dunit3.c_str(), com, &status);
    if (bunit!="") fits_update_key_str(fptr, "BUNIT", bunit.c_str(), com, &status);

    if (bmaj!=0) fits_update_key_dbl(fptr, "BMAJ", bmaj, 10, com, &status);
    if (bmin!=0) fits_update_key_dbl(fptr, "BMIN", bmin, 10, com, &status);
    fits_update_key_dbl(fptr, "BPA", bpa, 10, com, &status);
    if (btype!="NONE" && btype!="") fits_update_key_str(fptr, "BTYPE", btype.c_str(), com, &status);
    //fits_update_key_flt(fptr, "BZERO", bzero, 10, com, &status);
    //fits_update_key_flt(fptr, "BSCALE", bscale, 10, com, &status);
    //fits_update_key_flt(fptr, "BLANK", blank, 10, com, &status);

    if (object!="NONE" && object!="") fits_update_key_str(fptr, "OBJECT", object.c_str(), com, &status);
    if (epoch!=0) fits_update_key_flt(fptr, "EPOCH", epoch, 10, com, &status);
    if (telescope!="NONE" && telescope!="") fits_update_key_str(fptr, "TELESCOP", telescope.c_str(), com, &status);
    if (freq0!=0) fits_update_key_dbl(fptr, "FREQ0", freq0, 10, com, &status);
    if (datamax!=0) fits_update_key_dbl(fptr, "DATAMAX", datamax, 10, com, &status);
    if (datamin!=0) fits_update_key_dbl(fptr, "DATAMIN", datamin, 10, com, &status);

    if (fullHead) {
        char *Keys = new char[100];
        for (uint i=0; i<keys.size(); i++) {
            status=0;
            bool towrite = false;
            int hist = keys[i].find("HISTORY");

            if(hist>=0) towrite=true;
            else {
                int found = keys[i].find("=");
                if (found>=0) {
                    char keyname [] = "                                                                                  ";
                    strncpy(keyname, keys[i].c_str(),found);
                    char card[100];
                    if (fits_read_card(fptr, keyname, card, &status)) towrite=true;
                }
            }
            if (towrite) {
                status=0;
                strcpy(Keys, keys[i].c_str());
                fits_write_record(fptr, Keys, &status);
            }
        }
        delete [] Keys;
    }


    fits_report_error(stderr, status);

}


template <class T>
inline bool FitsHeader::read_keyword(std::string keyword, T &key, bool err) {

    int datatype=-1;

    if (std::is_same<T, int>::value) datatype=TINT;
    else if (std::is_same<T, long>::value) datatype=TLONG;
    else if (std::is_same<T, float>::value) datatype=TFLOAT;
    else if (std::is_same<T, double>::value) datatype=TDOUBLE;
    else if (std::is_same<T, long>::value) datatype=TSTRING;
    else {
        std::cout << "Error: unknown type of keyword "<< keyword << std::endl;
        return false;
    }

    fitsfile *fptr;
    int status=0;
    char comment[72];

    if (fits_open_file(&fptr, fitsname.c_str(), READONLY, &status)) {
        fits_report_error(stderr, status);
        return false;
    }
    status = 0;
    if (fits_read_key(fptr, datatype, keyword.c_str(), &key, comment, &status)) {
        if (err) fits_report_error(stderr, status);
        return false;
    }
    status = 0;
    if (fits_close_file(fptr, &status))
        fits_report_error(stderr, status);

    return true;
}


inline int FitsHeader::pixToWCS(const double *pix, double *world, size_t npts) {

    ///   Uses wcs to convert the three-dimensional pixel positions referenced 
    ///   by pix to world coordinates, which are placed in the array world[].
    ///   pix is assumed to hold the positions of npts points.
    ///   Offsets these pixel values by 1 to account for the C arrays being 
    ///   indexed to 0.
    /// 
    ///   \param pix The array of pixel coordinates.
    ///   \param world The returned array of world coordinates.
    ///   \param npts The number of distinct pixels in the arrays.

    int naxis=wcs->naxis,status;
    int specAxis = wcs->spec;
    if(specAxis<0) specAxis=2;
    if(specAxis>=naxis) specAxis = naxis-1;

    // correct from 0-indexed to 1-indexed pixel array
    // Add entries for any other axes that are present, 
    // keeping the order of pixel positions the same

    double *newpix = new double[naxis*npts];
    for(int pt=0;pt<npts;pt++){
      for(int i=0;i<naxis;i++) newpix[pt*naxis+i] = 1.;
      newpix[pt*naxis+wcs->lng]  += pix[pt*3+0];
      newpix[pt*naxis+wcs->lat]  += pix[pt*3+1];
      newpix[pt*naxis+specAxis] += pix[pt*3+2];
    }

    int    *stat      = new int[npts];
    double *imgcrd    = new double[naxis*npts];
    double *tempworld = new double[naxis*npts];
    double *phi       = new double[npts];
    double *theta     = new double[npts];
    status=wcsp2s(wcs, npts, naxis, newpix, imgcrd, phi, theta, tempworld, stat);
    
    if(status>0){
        std::cerr << "\nCannot convert to wcs. WCS error code = " << status
                  << ": stat="<<stat[0] << " : " << wcs_errmsg[status] << std::endl;
    }
    else {
        //return just the spatial/velocity information, keeping the
        //  order of the pixel positions the same.
        for(int pt=0;pt<npts;pt++){
            world[pt*3+0] = tempworld[pt*naxis + wcs->lng];
            world[pt*3+1] = tempworld[pt*naxis + wcs->lat];
            world[pt*3+2] = tempworld[pt*naxis + specAxis];
        }
    }
    
    delete [] stat;
    delete [] imgcrd;
    delete [] tempworld;
    delete [] phi;
    delete [] theta;
    delete [] newpix;
    return status;
    
}


int FitsHeader::wcsToPix(const double *world, double *pix, size_t npts) 
{
  ///  @details
  ///   Uses wcs to convert the three-dimensional world coordinate position 
  ///    referenced by world to pixel coordinates, which are placed in the
  ///    array pix[].
  ///   world is assumed to hold the positions of npts points.
  ///   Offsets the pixel values by 1 to account for the C arrays being 
  ///    indexed to 0.
  /// 
  /// \param wcs The wcsprm struct containing the WCS information.
  /// \param world The array of world coordinates.
  /// \param pix The returned array of pixel coordinates.
  /// \param npts The number of distinct pixels in the arrays.

  int naxis=wcs->naxis,status=0;
  int specAxis = wcs->spec;
  if(specAxis<0) specAxis=2;
  if(specAxis>=naxis) specAxis = naxis-1;

  // Test to see if there are other axes present, eg. stokes
  if(wcs->naxis>naxis) naxis = wcs->naxis;

  // Add entries for any other axes that are present, keeping the 
  //   order of pixel positions the same
  double *tempworld = new double[naxis*npts];
  for(int pt=0;pt<npts;pt++){
    for(int axis=0;axis<naxis;axis++) 
      tempworld[pt*naxis+axis] = wcs->crval[axis];
    tempworld[pt*naxis + wcs->lng]  = world[pt*3+0];
    tempworld[pt*naxis + wcs->lat]  = world[pt*3+1];
    tempworld[pt*naxis + specAxis] = world[pt*3+2];
  }

  int    *stat   = new int[npts];
  double *temppix = new double[naxis*npts];
  double *imgcrd = new double[naxis*npts];
  double *phi    = new double[npts];
  double *theta  = new double[npts];
  status=wcss2p(wcs,npts,naxis,tempworld,phi,theta,
        imgcrd,temppix,stat);
  if(status>0){
      std::cerr << "\nCannot convert to wcs. WCS error code = " << status
                << ": stat="<<stat[0] << " : " << wcs_errmsg[status] << std::endl;
  }
  else{
    // correct from 1-indexed to 0-indexed pixel array 
    //  and return just the spatial/velocity information, 
    //  keeping the order of the pixel positions the same.
    for(int pt=0;pt<npts;pt++){
      pix[pt*naxis + 0] = temppix[pt*naxis + wcs->lng] - 1.;
      pix[pt*naxis + 1] = temppix[pt*naxis + wcs->lat] - 1.;
      pix[pt*naxis + 2] = temppix[pt*naxis + specAxis] - 1.;
    }
  }

  delete [] stat;
  delete [] imgcrd;
  delete [] temppix;
  delete [] phi;
  delete [] theta;
  delete [] tempworld;
  return status;
}

#endif
