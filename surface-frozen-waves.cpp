#include <iostream>
#include <string>
#include <complex>
#include <cmath>
#include <iomanip>
#include <ctime>
#include <cstdlib>
#include "./../bessel-library/bessel-library.hpp" // A library to evaluate Bessel functions available
                                                  //  on https://github.com/jodesarro/bessel-library
#include "./../mtx-library/mtx-library.cpp" // A library to handle mtx files available
                                            //  on https://github.com/jodesarro/mtx-library

#ifndef M_PI
#define M_PI  3.1415926535897932384626433
#endif

#ifndef M_2PI
#define M_2PI 6.2831853071795864769252867
#endif

using namespace std;

struct SfwParameters
{
    complex<double> n;
    double l0, Q, L, H;
    int N, P;
    char* slice_axis_name = new char[1];
    double slice_axis_value, xmin, xmax, ymin, ymax, zmin, zmax;
    int xpoints, ypoints, zpoints, is_nonresistant;
    string file_in, file_out;
    long int start_time;
};

void log_information()
{
    puts("\n| SURFACE-FROZEN-WAVES | 21 SEP 2020 | @JODESARRO |\n");
    puts("DESCRIPTION");
    puts("Intensity calculation of ordinary Surface Frozen Waves (SFWs) "
        "with height and length proportional to F(x,z) and LFWs equally "
        "spaced with the same L and N. Documentation and latest release "
        "available on 'https://github.com/jodesarro/surface-frozen-waves'.\n");
}

void log_wrong_parameters_amount()
{
    puts("ERROR");
    puts("Error in processing parameters. Make sure to pass all "
    "parameters as argument.");
    puts("The program finished.");
}

void log_input_parameters( SfwParameters sfw )
{
    puts("INPUT PARAMETERS");
    cout << "nr = " << real(sfw.n) << endl;
    cout << "ni = " << imag(sfw.n) << endl;
    cout << "l0 = " << sfw.l0 << endl;
    cout << "Q = " << sfw.Q << endl;
    cout << "N = " << sfw.N << endl;
    cout << "L = " << sfw.L << endl;
    cout << "P = " << sfw.P << endl;
    cout << "slice_axis_name = " << sfw.slice_axis_name[0] << endl;
    cout << "slice_axis_value = " << sfw.slice_axis_value << endl;
    cout << "xmin = " << sfw.xmin << endl;
    cout << "xmax = " << sfw.xmax << endl;
    cout << "xpoints = " << sfw.xpoints << endl;
    cout << "ymin = " << sfw.ymin << endl;
    cout << "ymax = " << sfw.ymax << endl;
    cout << "ypoints = " << sfw.ypoints << endl;
    cout << "zmin = " << sfw.zmin << endl;
    cout << "zmax = " << sfw.zmax << endl;
    cout << "zpoints = " << sfw.zpoints << endl;
    cout << "file_in = " << sfw.file_in << endl;
    cout << "file_out = " << sfw.file_out << endl << endl;
}

void log_progress_calculating()
{
    puts("PROGRESS");
    puts("Calculating...");
}

int main(int argc, char *argv[])
{

    log_information();

    SfwParameters sfw;
    if ( argc != 22 )
    {
        log_wrong_parameters_amount();
        exit(1);
    }
    else
    {
        sfw.start_time = time(NULL);
        sfw.n = complex<double> ( atof(argv[1]), atof(argv[2]) );
        sfw.l0 = atof(argv[3]);
        sfw.Q  = atof(argv[4]);
        sfw.N  = atoi(argv[5]);
        sfw.L  = atof(argv[6]);
        sfw.P  = atoi(argv[7]);
        sfw.slice_axis_name = argv[8];
        sfw.slice_axis_value = atof(argv[9]);
        sfw.xmin = atof(argv[10]);
        sfw.xmax = atof(argv[11]);
        sfw.xpoints = atoi(argv[12]);
        sfw.ymin = atof(argv[13]);
        sfw.ymax = atof(argv[14]);
        sfw.ypoints = atoi(argv[15]);
        sfw.zmin = atof(argv[16]);
        sfw.zmax = atof(argv[17]);
        sfw.zpoints = atoi(argv[18]);
        sfw.is_nonresistant = atoi(argv[19]);
        sfw.file_in  = argv[20];
        sfw.file_out = argv[21];
    }

    log_input_parameters( sfw );
 
    log_progress_calculating();

    // === READ MTX FILE =============================================
    int pxH, pxL, pxHL;
    mtx_get_dimensions( sfw.file_in, pxH, pxL );
    pxHL = pxH*pxL;
    double * F = new double [ pxHL ];
    mtx_read_to_array( sfw.file_in, F );

    // === WAVE VECTORS ==============================================
    int N2M1 = sfw.N*2+1;
    complex<double> kk = (M_2PI*sfw.n/sfw.l0)*(M_2PI*sfw.n/sfw.l0);
    complex<double> * b = new complex<double> [ N2M1 ]; // beta
    complex<double> * h = new complex<double> [ N2M1 ]; // krho
    // In what follows, iq=q+N, q=iq-N, [iq]->q 
    for ( int iq=0; iq<N2M1; iq++ )
    {
        b[iq] = complex<double> ( sfw.Q + M_2PI*(iq-sfw.N)/sfw.L, ( sfw.Q + M_2PI*(iq-sfw.N)/sfw.L )*imag(sfw.n)/real(sfw.n) );
        h[iq] = sqrt( kk - b[iq]*b[iq] );
    }

    // === COEFFICIENTS Aqps =========================================
    complex<double> * Aqp;
    if ( sfw.P < pxH )
    {
        complex<double> n_temp = ( sfw.is_nonresistant == 1 ) ? complex<double>( real(sfw.n), 0 ) : sfw.n; 
        Aqp = new complex<double> [ sfw.P*N2M1 ];

        int pxLm1 = pxL-1;

        for ( int ip=0; ip<sfw.P; ip++ )
        {
            // Variable that relates pixels in x of input image with number P of LFW
            int pxHip_P = int ( floor(pxH*ip/sfw.P) );

            for ( int iq=0; iq<N2M1; iq++ )
            {
                // The following is iteration in j that represents pixels in z. It is an
                // approximation of the Aqp integral by trapezoidal method for equally spaced z.
                // Each j is associated to a point in z (or a pixel in zj of input image).
                complex<double> temp = 0.5*( F[ pxHip_P + pxH*pxLm1 ]*exp( sfw.L*sfw.Q*imag(n_temp)/real(n_temp) ) + F[ pxHip_P + 0 ] );
                for ( int j=1; j<pxLm1; j++ )
                {
                    temp += F[ pxHip_P + pxH*j ]*exp( complex<double>(0., -M_2PI*(iq-sfw.N)*j/pxLm1) )
                            * exp( (sfw.Q*imag(n_temp)/real(n_temp))*sfw.L*j/pxLm1 );                    
                }
                Aqp[ iq + N2M1*ip ] = temp/( double (pxLm1) );
            }
        }
    }
    else // It is an optmization for pxH < P, Aqp -> A[q][pxh]
    {
        complex<double> n_temp = ( sfw.is_nonresistant == 1 ) ? complex<double>( real(sfw.n), 0 ) : sfw.n; 
        Aqp = new complex<double> [ pxH*N2M1 ];

        int pxLm1 = pxL-1;

        for ( int ipxh=0; ipxh<pxH; ipxh++ )
        {

            for ( int iq=0; iq<N2M1; iq++ )
            {
                // The following is iteration in j that represents pixels in z. It is an
                // approximation of the Aqp integral by trapezoidal method for equally spaced z.
                // Each j is associated to a point in z (or a pixel in zj of input image).
                complex<double> temp = 0.5*( F[ ipxh + pxH*pxLm1 ]*exp( sfw.L*sfw.Q*imag(n_temp)/real(n_temp) ) + F[ ipxh + 0 ] );
                for ( int j=1; j<pxLm1; j++ )
                {
                    temp += F[ ipxh + pxH*j ]*exp( complex<double>(0., -M_2PI*(iq-sfw.N)*j/pxLm1) )
                            * exp( (sfw.Q*imag(n_temp)/real(n_temp))*sfw.L*j/pxLm1 );                    
                }
                Aqp[ iq + N2M1*ipxh ] = temp/( double (pxLm1) );
            }
        }
    }
    delete[] F;

    // === DISPLACEMENT IN X OF EACH LFW RELATIVE TO THE ORIGIN ======
    double * x0 = new double [ sfw.P ];
    double * x0x0 = new double [ sfw.P ];
    sfw.H = sfw.L * double(pxH)/double(pxL); 
    // In what follows, ip=p-1, p=ip+1, [ip]->p 
    for ( int ip=0; ip<sfw.P; ip++ )
    {
        // Linear and homogeneous distribution in x according to P and H.
        x0[ip] = ( double(ip) )*sfw.H/( double(sfw.P) );
        x0x0[ip] = x0[ip]*x0[ip];
    }

    // === INTENSITY OF SFW ==========================================
    double * sfwi;
    int mtx_rows, mtx_columns;
    if ( sfw.slice_axis_name[0] == 'y' || sfw.slice_axis_name[0] == 'Y' )
    {     
        sfw.ymin = 0;
        sfw.ymax = 0;
        sfw.ypoints = 0;
        sfwi = new double [ sfw.xpoints*sfw.zpoints ];
        double dx = (sfw.xmax-sfw.xmin)/( double(sfw.xpoints) );
        double dz = (sfw.zmax-sfw.zmin)/( double(sfw.zpoints) );
        double var_yy = sfw.slice_axis_value*sfw.slice_axis_value;

        // Longitudinal iteration
        for ( int iz=0; iz<sfw.zpoints; iz++ )
        {
            complex<double> var_iz = complex<double> (0., sfw.zmin + ( double(iz) ) * dz);

            // Radial iteration
            for ( int ix=0; ix<sfw.xpoints; ix++ )
            {
                double var_x = sfw.xmin + ( double(ix) ) * dx;
                double var_rhorho = var_x*var_x + var_yy;
                double var_2rho_cosphi = 2.0*var_x;
                complex<double> temp = complex<double> (0.0, 0.0);

                for ( int ip=0; ip<sfw.P; ip++ )
                {
                    if ( sfw.P < pxH ) // pxH >= P, Aqp -> A[q][p]
                    {
                        for ( int iq=0; iq<N2M1; iq++ )
                        {
                            temp += Aqp[ iq + N2M1*ip ]
                                    * bessel::cyl_j0( h[iq]*sqrt(var_rhorho+x0x0[ip]-var_2rho_cosphi*x0[ip]) )
                                    * exp( b[iq]*var_iz );
                        }
                    }
                    else // pxH < P, Aqp -> A[q][pxh]
                    {
                        // Variable that relates pixels in x of input image with number P of LFW
                        int pxHip_P = int ( floor(pxH*ip/sfw.P) );

                        for ( int iq=0; iq<N2M1; iq++ )
                        {
                            temp += Aqp[ iq + N2M1*pxHip_P ]
                                    * bessel::cyl_j0( h[iq]*sqrt(var_rhorho+x0x0[ip]-var_2rho_cosphi*x0[ip]) )
                                    * exp( b[iq]*var_iz );
                        }
                    }
                }
                sfwi[ix + sfw.xpoints*iz] = abs(temp)*abs(temp);
            }
        }
        mtx_rows = sfw.xpoints;
        mtx_columns = sfw.zpoints;
    }
    else if ( sfw.slice_axis_name[0] == 'z' || sfw.slice_axis_name[0] == 'Z' )
    {
        sfw.zmin = 0;
        sfw.zmax = 0;
        sfw.zpoints = 0;
        sfwi = new double [ sfw.xpoints*sfw.ypoints ];
        double dx = (sfw.xmax-sfw.xmin)/( double(sfw.xpoints) );
        double dy = (sfw.ymax-sfw.ymin)/( double(sfw.ypoints) );
        complex<double> var_iz = complex<double> (0., sfw.slice_axis_value);

        // Transversal iteration
        for ( int iy=0; iy<sfw.ypoints; iy++ )
        {
            double var_yy = (sfw.ymin + ( double(iy) ) * dy)
                          * (sfw.ymin + ( double(iy) ) * dy);

            // Radial iteration 
            for ( int ix=0; ix<sfw.xpoints; ix++ )
            {
                double var_x = sfw.xmin + ( double(ix) ) * dx;
                double var_rhorho = var_x*var_x + var_yy;
                double var_2rho_cosphi = 2.0*var_x;
                complex<double> temp = complex<double> (0.0, 0.0);

                for ( int ip=0; ip<sfw.P; ip++ )
                {
                    if ( sfw.P < pxH ) // pxH >= P, Aqp -> A[q][p]
                    {
                        for ( int iq=0; iq<N2M1; iq++ )
                        {
                            temp += Aqp[ iq + N2M1*ip ]
                                    * bessel::cyl_j0( h[iq]*sqrt(var_rhorho+x0x0[ip]-var_2rho_cosphi*x0[ip]) )
                                    * exp( b[iq]*var_iz );
                        }
                    }
                    else // pxH < P, Aqp -> A[q][pxh]
                    {
                        // Variable that relates pixels in x of input image with number P of LFW
                        int pxHip_P = int ( floor(pxH*ip/sfw.P) );
                       
                        for ( int iq=0; iq<N2M1; iq++ )
                        {
                            temp += Aqp[ iq + N2M1*pxHip_P ]
                                    * bessel::cyl_j0( h[iq]*sqrt(var_rhorho+x0x0[ip]-var_2rho_cosphi*x0[ip]) )
                                    * exp( b[iq]*var_iz );
                        }
                    }
                }
                sfwi[ix + sfw.xpoints*iy] = abs(temp)*abs(temp);
            }
        }
        mtx_rows = sfw.xpoints;
        mtx_columns = sfw.ypoints;
    }
    else if ( sfw.slice_axis_name[0] == 'x' || sfw.slice_axis_name[0] == 'X' )
    {
        sfw.xmin = 0;
        sfw.xmax = 0;
        sfw.xpoints = 0;
        sfwi = new double [ sfw.ypoints*sfw.zpoints ];
        double dy = (sfw.ymax-sfw.ymin)/( double(sfw.ypoints) );
        double dz = (sfw.zmax-sfw.zmin)/( double(sfw.zpoints) );
        double var_xx = sfw.slice_axis_value*sfw.slice_axis_value;
        double var_2rho_cosphi = 2.0*sfw.slice_axis_value;

        // Longitudinal iteration
        for ( int iz=0; iz<sfw.zpoints; iz++ )
        {
            complex<double> var_iz = complex<double> (0., sfw.zmin + ( double(iz) ) * dz);

            // Transversal iteration
            for ( int iy=0; iy<sfw.ypoints; iy++ )
            {
                double var_rhorho = var_xx
                 + (sfw.ymin + ( double(iy) ) * dy)*(sfw.ymin + ( double(iy) ) * dy); // rho²=x²+y²
                complex<double> temp = complex<double> (0.0, 0.0);

                for ( int ip=0; ip<sfw.P; ip++ )
                {
                    if ( sfw.P < pxH ) // pxH >= P, Aqp -> A[q][p]
                    {
                        for ( int iq=0; iq<N2M1; iq++ )
                        {
                            temp += Aqp[ iq + N2M1*ip ]
                                    * bessel::cyl_j0( h[iq]*sqrt(var_rhorho+x0x0[ip]-var_2rho_cosphi*x0[ip]) )
                                    * exp( b[iq]*var_iz );
                        }
                    }
                    else // pxH < P, Aqp -> A[q][pxh]
                    {
                        // Variable that relates pixels in x of input image with number P of LFW
                        int pxHip_P = int ( floor(pxH*ip/sfw.P) );
                        
                        for ( int iq=0; iq<N2M1; iq++ )
                        {
                            temp += Aqp[ iq + N2M1*pxHip_P ]
                                    * bessel::cyl_j0( h[iq]*sqrt(var_rhorho+x0x0[ip]-var_2rho_cosphi*x0[ip]) )
                                    * exp( b[iq]*var_iz );
                        }
                    }
                }
                sfwi[iy + sfw.ypoints*iz] = abs(temp)*abs(temp);
            }
        }
        mtx_rows = sfw.ypoints;
        mtx_columns = sfw.zpoints;
    }
    
    delete[] Aqp;
    delete[] x0;
    delete[] x0x0;
    delete[] b;
    delete[] h;

    // === WRITE MTX FILE ============================================
    string mtx_comments;
    mtx_comments = "nr ni l0 Q N L P H slice_axis_name slice_axis_value xmin xmax xpoints ymin ymax ypoints zmin zmax zpoints is_nonresistant elapsed_time\n" "%" + to_string(real(sfw.n)) + " " + to_string(imag(sfw.n)) + " " + to_string(sfw.l0) + " " + to_string(sfw.Q) + " " + to_string(sfw.N) + " " + to_string(sfw.L) + " " + to_string(sfw.P) + " " + to_string(sfw.H) + " " + sfw.slice_axis_name[0] + " " + to_string(sfw.slice_axis_value) + " " + to_string(sfw.xmin) + " " + to_string(sfw.xmax) + " " + to_string(sfw.xpoints) + " " + to_string(sfw.ymin) + " " + to_string(sfw.ymax) + " " + to_string(sfw.ypoints) + " " + to_string(sfw.zmin) + " " + to_string(sfw.zmax) + " " + to_string(sfw.zpoints) + " " + to_string(sfw.is_nonresistant) + " " + to_string( time(NULL) - sfw.start_time )  + "\n";
    mtx_write_from_array( sfw.file_out, mtx_comments, mtx_rows, mtx_columns, sfwi );

    delete[] sfwi;
    
    cout << "Done.\nElapsed time: " << ( time(NULL) - sfw.start_time ) << " seconds." << endl;
}