#include <iostream>
#include <fitsarray.hh>
#include <pvslice.hh>


int main (int argc, char *argv[]) {
    
    if (argc!=6 && argc!=7) {
        std::cout << "PVCREATOR usage: \n\n"
                  << "   1) pvcreator <inp_cube> <out_pv> <x0> <y0> <angle> \n"
                  << "   2) pvcreator <inp_cube> <out_pv> <x1> <y1> <x2> <y2> \n\n"
                  << "In 1) the slice is defined by the point (x0,y0) and angle (N->W).\n"
                  << "In 2) the slice is defined by two points (x1,y1) and (x2,y2)\n";
        return EXIT_FAILURE;
    }
    else {
        std::string infile = argv[1];
        std::string outfile = argv[2];
        FitsCube<float> c(infile);
        PvSlice<float> *s;
        
        if (argc==6) {
            float x0 = atof(argv[3]), y0 = atof(argv[4]), a = atof(argv[5]);
            s =  new PvSlice<float>(x0,y0,a);  
        }
        else if (argc==7) {
            float x1 = atof(argv[3]), y1 = atof(argv[4]), x2 = atof(argv[5]), y2 = atof(argv[6]);
            s =  new PvSlice<float>(x1,y1,x1,y2);  
        }
        
        std::cout << "Extracting PV and writing it to " << outfile << " ..." << std::flush;
        s->slice(&c);
        s->fitswrite(outfile,true);
        std::cout << "Done! " << std::endl;
    
    }
    
    return EXIT_SUCCESS;
}