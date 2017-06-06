//------------------------------------------------------------------
// pvslicer.hh
// Definition of PvSlice class
//------------------------------------------------------------------  

#ifndef PVSLICE_HH
#define PVSLICE_HH

#include <iostream>
#include <array.hh>
#include <fitsarray.hh>
#include <cstring>
#include <cmath>

//////////////////////////////////////////////////////////////////////////
// PvSlicer class to extract a position-velocity slice from a datacube
// The slice is always a straigh line and can be defined:
//  1) Through two points (constructor #1)
//  2) Through one point and an angle (constructor #2)
//
// In the first case, the slice can include just a part of the datacube,
// while in the second the entire cube is sliced.
// Usage: call one of the constructors and the slice(FitsCube).
//////////////////////////////////////////////////////////////////////////
template <class T>
class PvSlice : public FitsImage<T>
{
public:
    PvSlice(float x1, float y1, float x2, float y2) : x1(x1), x2(x2), y1(y1), y2(y2) {}
    PvSlice(float x0, float y0, float angle) : x0(x0), y0(y0), angle(angle), isAngle(true) {}
    ~PvSlice(){if (x_locus) delete [] x_locus; if (y_locus) delete [] y_locus;}
    bool slice (FitsCube<T> *c);
    
private:
    size_t xpix, ypix, zpix;        //< Size of the datacube to slice
    float  x1, x2, y1, y2;          //< Defining slice with two points
    float  x0, y0, angle;           //< Defining slice with a point and a angle (from North->West)
    bool   isAngle;                 //< Whether slice is defined with angle
    float  *x_locus, *y_locus;      //< The slice
    size_t num_points;              //< Number of pixels along the slice
    
    
    float weight (float x, float y, float cx, float cy) {return fabs((1-(x-cx))*(1-(y-cy)));}  
    bool  define_slice(int x1, int y1, int x2, int y2);
    bool  check_bounds (int *blx, int *bly, int *Trx, int *Try);
    bool  pvslice (FitsCube<T> *c);
    void  define_header(FitsCube<T> *c);
};


template <class T>
bool PvSlice<T>::slice(FitsCube<T> *c) {
    
    // Front-end function to slice the cube and extract the position-velocity.
    // Slice is made through the data provided in the constructor.
    
    xpix = c->DimX();
    ypix = c->DimY();
    zpix = c->DimZ();
    
    if (!(xpix*ypix*zpix > 0)) {
        fprintf (stderr, "PvSlice ERROR: input cube dimensions are wrong.\n");
        return false;
    }
    
    // If the slice is defined through a point and a angle, the slice is always
    // taken across the entire cube, so calculate intercepts at 0 and xsize-1
    if (isAngle) {
        if (x0<0 || x0>=xpix || y0<0 || y0>=ypix) {
            fprintf (stderr, "PvSlice ERROR: center (x0,y0) is outside the cube.\n");
            return false;
        }
        while (angle>=360) angle-=360;
        while (angle<=-360) angle+=360;
        if (angle==0 || angle==180) {
            x1 = x2 = x0;
            y1 = 0; y2 = ypix-1;
        }
        else {
            float theta = (angle+90)*M_PI/180.;
            x1 = 0;      y1 = tan(theta)*(x1-x0)+y0;
            x2 = xpix-1; y2 = tan(theta)*(x2-x0)+y0;
        }
    }
        
    // Get the slice locus 
    if(!define_slice(x1,y1,x2,y2)) return false;
    
    if (num_points>0) {
        // Setting the PV array
        size_t dimen[2] = {num_points,zpix};
        this->setFitsArray(dimen);
        
        // Extract slice from cube 
        if (!pvslice(c)) return false;
        define_header(c);
    }
    
    return true;
}


template <class T>
bool PvSlice<T>::define_slice(int x1, int y1, int x2, int y2) {

    // Determine a locus of pixels in a slice line, from two endpoints 
    
    int    blx, bly, Trx, Try;
    float  theta, ctheta, stheta;

    if (!check_bounds(&blx,&bly,&Trx,&Try) ) return false;

    /* now calculate slice length */
    float dx    = (float) (Trx-blx+1);
    float dy    = (float) (Try-bly+1);
    if (blx==Trx) {
        theta = M_PI/2.0; ctheta=0.0; stheta=1.0;
        if (bly>Try) {theta = 3.0*M_PI/2.0; stheta=-1.0;}
    }
    else {theta = atan2(dy,dx); ctheta = cos(theta); stheta = sin(theta);}

    float slicelength = sqrt(dx*dx+dy*dy);
    num_points = lround(slicelength);

    if (x_locus) delete [] x_locus;
    if (y_locus) delete [] y_locus;
    x_locus = new float[num_points];
    y_locus = new float[num_points];
    
    for (int i=0; i<num_points; i++) {
        x_locus[i] = i * ctheta + blx;
        y_locus[i] = i * stheta + bly;
    }

    return true;
} 


template <class T>
bool PvSlice<T>::check_bounds (int *blx, int *bly, int *Trx, int *Try) {
    
    // Checks bounds & if endpoints are outside the area of cube's front face,
    // computes overlapping section of slice 

    float f, x, y;
    int swapped=0, v[5], w[5];
    float ix1 = x1, iy1 = y1,  ix2 = x2, iy2 = y2;

    /* do everything assuming first point is left of second, swap back later */
    if (ix2 < ix1) {
        ix1 = x2; ix2 = x1; iy1 = y2; iy2 = y1; swapped=1;
    } 
    else {
        if (ix1 == ix2 && iy1 > iy2) {
            iy1 = y2; iy2 = y1; swapped=2;
        }
    }
   
    *blx = ix1; *bly = iy1; *Trx = ix2; *Try = iy2;

    // Handle no-overlap cases: the extended slice line may overlap eventually,
    // but the piece within the endpoints does not. So stop here.
    if ( (ix1 < 0 && ix2 < 0) || (ix1 >= xpix && ix2 >= xpix) ) {
        fprintf(stderr, "PvSlice ERROR: Slice does not intercept cube.\n"); return false;
    }
    if ( (iy1 < 0 && iy2 < 0) || (iy1 >= ypix && iy2 >= ypix) ) {
        fprintf(stderr, "PvSlice ERROR: Slice does not intercept cube.\n"); return false;
    }
    // Single pixel case
    if ( ix1 == ix2 && iy1 == iy2 ) {
        fprintf(stderr, "PvSlice ERROR: Line segment too short\n"); return false;
    }

    // Handle the case of a vertical slice 
    if ( ix1 == ix2 ) {

        if (ix1 < 0 || ix1 >= xpix) {
            fprintf(stderr, "PvSlice ERROR: Slice does not intercept cube.\n"); return false;
        }
        else {
            *blx = ix1; *Trx = ix2;
            // Now check if either of the y coords are good 
            if (iy1 < 0 || iy1 >= ypix) { *bly = 0;}
            if (iy2 < 0 || iy2 >= ypix) { *Try = ypix-1;}
        }
    }
    else {
        // Handle the case of a horizontal slice
        if ( iy1 == iy2 ) {

            if (iy1 < 0 || iy1 >= ypix) {
                fprintf(stderr, "PvSlice ERROR: Slice does not intercept cube.\n"); return false;
            }
            else {
                *bly = iy1; *Try = iy2;
                // Now check if either of the x coords are good 
                if (ix1 < 0 || ix1 >= xpix) { *blx = 0;}
                if (ix2 < 0 || ix2 >= xpix) { *Trx = xpix-1;}
            }
        }
        else {
            // Calculate intercepts with x=0, y=0, x=xpix-1, y=ypix-1
            // store x-intercepts in v[] and y-intercepts in w[]       
            f = ( (float)(iy2) - (float)(iy1) ) / ( (float)(ix2) - (float)(ix1) );

            v[1]=0; w[2]=0; v[3]=xpix-1; w[4]=ypix-1;
            x=0.0;    y=(float)(iy1) + (x-(float)(ix1))*f;          w[1]=(int)(y+0.5);
            y=0.0;    x=(float)(ix1) + (y-(float)(iy1))/f;          v[2]=(int)(x+0.5);
            x=(float)(xpix-1); y=(float)(iy1) + (x-(float)(ix1))*f; w[3]=(int)(y+0.5);
            y=(float)(ypix-1); x=(float)(ix1) + (y-(float)(iy1))/f; v[4]=(int)(x+0.5);


            // Things are different depending on if slope is +ve or -ve 
            // For either, there are six cases of where the intercepts of the
            // slice line (extended to infinity in both directions) lie w.r.t. the
            // four lines noted above
            if ( f >= 0.0 ) {
                if (v[2] < 0 ) {
                    if (w[1] < ypix) {
                        *blx = 0; *bly = w[1];
                        if ( w[3] < ypix ) {*Trx = xpix-1; *Try = w[3];}
                        else{*Trx = v[4]; *Try = ypix-1;}
                    }
                    else {fprintf(stderr, "PvSlice ERROR: Out of bounds\n"); return false;}
                }
                else{
                    if (v[2] < xpix) {
                        *bly = 0; *blx = v[2];
                        if (v[4] < xpix) {*Try = ypix-1; *Trx = v[4];}
                        else {*Try = w[3]; *Trx = xpix-1;}
                    }
                    else{fprintf(stderr, "PvSlice ERROR: Out of bounds\n"); return false;}
                }
            }
            else {
                if (v[4] < 0) {
                    if ( w[1] > 0 ) {
                        *blx = 0; *bly = w[1];
                        if (w[3] > 0 ) {*Trx = xpix-1; *Try = w[3];}
                        else {*Trx = v[2]; *Try = 0;}
                    }
                    else {fprintf(stderr, "PvSlice ERROR: Out of bounds\n"); return false;}
                }
                else {
                    if (v[4] < xpix) {
                        *blx = v[4]; *bly = ypix-1;
                        if ( v[2] < xpix ) {*Trx = v[2]; *Try = 0;}
                        else {*Trx = xpix-1; *Try = w[3];}
                    }
                    else {fprintf(stderr, "PvSlice ERROR: Out of bounds\n"); return false;}
                }
            } 
        } 
    } 
    
    // now check to see if either endpoint NEEDS to be replaced:
    // if they are legal return the input value 
    if ( (ix1 >= 0 && ix1 < xpix) && (iy1 >= 0 && iy1 < ypix) ) {
        *blx = ix1; *bly = iy1;
    }
    if ( (ix2 >= 0 && ix2 < xpix) && (iy2 >= 0 && iy2 < ypix) ) {
        *Trx = ix2; *Try = iy2;
    }

    // Swap ends back if necessary
    if (swapped > 0) {
        ix1 = *blx; *blx = *Trx; *Trx = ix1;
        iy1 = *bly; *bly = *Try; *Try = iy1;
    }

    return true;

} 


template <class T>
bool PvSlice<T>::pvslice (FitsCube<T> *cube) {

    // Extract pixels from a cube, along an input locus and write the PV in 
    // the main array. Result is antialiased.
    // 
    // The cube is assumed to have axes in x,y,v order.
    // The output array has the same spatial scale as input cube, ie
    // the scale is assumed to be the same for both spatial axes,
    // and velocity pixels are given the same width as the channel spacing.
    // Output is a weighted sum over all pixels nearby the point where the
    // slice locus passes, to reduce aliasing effects.
    //
    // There are four cases to consider for the neighbour pixels, depending
    // where within a pixel the hit occurs:
    //
    //         |-----------|-----------|     |-----------|-----------|
    //         |           |           |	   |           |           |
    //    2.0  -     3     |     2     | 	   |     3     |     2     |
    //         |           |           |	   |       xy  |           |
    //         |-----------|-----------|	   |-----------|-----------|
    //         |           | xy        |	   |           |           |
    //    1.0  -     4     |     1     |	   |     4     |     1     |
    //         |           |           |	   |           |           |
    //         |-----|-----|-----|-----|	   |-----------|-----------|
    //              1.0         2.0
    //
    //         |-----------|-----------|     |-----------|-----------|
    //         |           |           |	   |           |           |
    //         |     3     |     2     | 	   |     3     |     2     |
    //         |           | xy        |	   |           |           |
    //         |-----------|-----------|	   |-----------|-----------|
    //         |           |           |	   |       xy  |           |
    //         |     4     |     1     |	   |     4     |     1     |
    //         |           |           |	   |           |           |
    //         |-----------|-----------|	   |-----------|-----------|
    //
    // All these can be handled by int(x+/-0.5), int(y+/-0.5)
    // Pixel coordinates are associated with centres of pixels. In pixels that 
    // don't have 4 neighbours, missing neighbours are aasigned a zero weight.

    
    int    xp, yp, xn, yn;
    float  xc, yc, xx, yy, sw, f;
    short MAXNB = 5;    // Number of neighbours for antialiasing
    

    FitsCube<float> c = *cube;
    if (num_points < 2) return false;
    
    float2D wt (MAXNB,num_points);
    int3D nb(2,MAXNB,num_points);

    // Find the neighbour pixels & antialiasing weights.
    // this is done for one channel only, then list applied to all channels
    for (int i = 0; i < num_points; i++) {
        xc = x_locus[i];       yc = y_locus[i];
        xp = (int)( xc+0.5 );  yp = (int)( yc+0.5 );
        // here I only do the 4 nearest pixels
        nb(0,0,i) = xp; nb(1,0,i) = yp;
        nb(0,1,i) = (int)(xc+0.5); nb(1,1,i) = (int)(yc-0.5); // neigbours
        nb(0,2,i) = (int)(xc+0.5); nb(1,2,i) = (int)(yc+0.5);
        nb(0,3,i) = (int)(xc-0.5); nb(1,3,i) = (int)(yc+0.5);
        nb(0,4,i) = (int)(xc-0.5); nb(1,4,i) = (int)(yc-0.5);

        // calculate the weight for each neighbour
        for (int j = 1; j < MAXNB; j++) {
            xp = nb(0,j,i); yp = nb(1,j,i);
            xx = (float) xp;   yy = (float) yp;
            if ( (xp >= 0) && (xp < xpix) && (yp >= 0) && (yp < ypix) )
                wt(j, i) = weight (xx, yy, xc, yc);
            else wt(j, i) = 0.0;
        }   
    }   

    
    for (int z = 0; z < zpix; ++z) {
        for (int i=0; i<num_points; i++) { /* start slice loop */
            f = 0.0; sw = 0.0;
            xp = nb(0,0,i);
            yp = nb(1,0,i);

            for (int j=1; j<MAXNB; j++) {
                xn = nb(0,j,i);
                yn = nb(1,j,i);
                /* The first test protects the second; it is possible for xn and yn
                to go out of range (neighbours of a pixel at edge of cube) */
                if ( wt(j,i) >0.0 && c(xn,yn,z)!=1.0e30) {
                    sw += wt(j,i);
                    f  += c(xn,yn,z) * wt(j,i);
                }
            } 

            if (sw>0.0) this->array(i,z) = f/sw;
            else {
                if ( xp >= 0 && xp<xpix && yp>=0 && yp<ypix ) this->array(i,z) = c(xn,yn,z);
                else this->array(i,z) = 0.0;
            }
        } 
    } 
    
    return true;
} 


template <class T>
void PvSlice<T>::define_header (FitsCube<T> *c) {
    
    // Defining the Header for the PV slice.
    
    ///* Calculating x-axis WCS
    double *coord_start = c->getXYphys(x_locus[0],y_locus[0]);
    double *coord_end   = c->getXYphys(x_locus[num_points-1],y_locus[num_points-1]);
    double coord_center[2] = {0.5*(coord_start[0]+coord_end[0]),0.5*(coord_start[1] + coord_end[1])};
    
    double ra_off  = (coord_start[0]-coord_center[0])*cos(coord_center[1]*M_PI/180.);
    double dec_off = coord_end[1] - coord_center[1];
    double position_angle = atan2(-ra_off, dec_off)*180./M_PI - 90.0;
    if (position_angle < 0) position_angle += 360.0;
    double first_coord  = -sqrt( ra_off*ra_off + dec_off*dec_off );
    
    double cdelt1 = fabs(2*first_coord/(num_points-1));
    double pcent = 0.5*num_points;
    
    if (isAngle) {
        // If the center is given, let's center the WCS on it.
        double dx_0 = sqrt( (x0-x_locus[0])*(x0-x_locus[0]) + (y0-y_locus[0])*(y0-y_locus[0]) ) ;
        double dx_t = sqrt( (x_locus[0]-x_locus[num_points-1])*(x_locus[0]-x_locus[num_points-1]) 
                          + (y_locus[0]-y_locus[num_points-1])*(y_locus[0]-y_locus[num_points-1]) );
        pcent =  dx_0/dx_t*num_points;
    }
    
    // Setting spatial coordinates
    this->Head().setCtype(0,"OFFSET");
    this->Head().setCrpix(0,pcent);
    this->Head().setCrval(0,0);
    this->Head().setCdelt(0,cdelt1);
    this->Head().setCunit(0,c->Head().CUNIT(0));
    
    // Setting spectral coordinates
    this->Head().setCtype(1,c->Head().CTYPE(2));    
    this->Head().setCrpix(1,c->Head().CRPIX(2));    
    this->Head().setCrval(1,c->Head().CRVAL(2));    
    this->Head().setCdelt(1,c->Head().CDELT(2));    
    this->Head().setCunit(1,c->Head().CUNIT(2));    
    
    // Setting other properties
    this->Head().setBmaj(c->Head().BMAJ());    
    this->Head().setBmin(c->Head().BMIN());    
    this->Head().setBpa(c->Head().BPA());    
    this->Head().setBunit(c->Head().BUNIT());
    this->Head().setBtype(c->Head().BTYPE());        
    this->Head().setEpoch(c->Head().EPOCH());    
    this->Head().setName(c->Head().OBJECT());    
    this->Head().setTelesc(c->Head().TELESCOPE());    
    
    std::string coordstr = "RA "+std::to_string(coord_center[0])+" deg, DEC "+std::to_string(coord_center[1])+" deg";
    this->Head().Keys().push_back("HISTORY BBAROLO PVSLICE: Center coordinates "+coordstr);
    this->setHeadDefined(true);
    
}


#endif
