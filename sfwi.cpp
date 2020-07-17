#include <iostream>
#include <fstream>
#include <string>
#include <complex>
#include <cmath>
#include <iomanip>
#include <ctime>
#include <cstdlib>
#include "./../bessel/bessel.hpp" // A library to evaluate Bessel functions
                                  //  of complex arguments available on
                                  //  https://github.com/jodesarro/bessel

#ifndef M_PI
#define M_PI 3.1415926535897932384626433
#endif

using namespace std;

void info();

int main(int argc, char *argv[])
{


    // === COUT INFO =================================================
    puts("\n| SFWI | 31 MAY 2020 | @JODESARRO |\n");
    puts("DESCRIPTION");
    puts("Intensity calculation of ordinary Surface Frozen Wave (SFW) "
        "with height and length proportional to F(x,z) and LFWs equally "
        "spaced with the same L and N. Documentation and latest release "
        "available on 'https://github.com/jodesarro/sfwi'.\n");
    // ===============================================================




    // === GENERAL ===================================================
    const double M_2PI = 2.0*M_PI;
    long int tempo = time(NULL);    // To count elapsed time          
    // ===============================================================




    // === COLLECTING THE INPUT PARAMETERS ===========================
    if ( argc != 22 )
    {
        puts("Error in processing parameters. Make sure to pass all "
        "parameters as argument.");
        puts("The program finished.");
        exit(1);
    }

    complex<double> n = complex<double>
     ( atof(argv[1]), atof(argv[2]) );
    double l0 = atof(argv[3]);
    double Q  = atof(argv[4]);
    double N  = atoi(argv[5]);
    double L  = atof(argv[6]);
    int P  = atoi(argv[7]);
    char * cut_axis = argv[8];         // Axis where the cut will be perpendicular to it
    double cut_value = atof(argv[9]);  // Value to cut the cut_axis
    double xmin = atof(argv[10]);
    double xmax = atof(argv[11]);
    int xpoints = atoi(argv[12]);
    double ymin = atof(argv[13]);
    double ymax = atof(argv[14]);
    int ypoints = atoi(argv[15]);
    double zmin = atof(argv[16]);
    double zmax = atof(argv[17]);
    int zpoints = atoi(argv[18]);
    int is_nonresistant = atoi(argv[19]);
    string file_in  = argv[20];
    string file_out = argv[21];
    // ===============================================================




    // === COUT PARAMETERS ===========================================
    puts("INPUT");
    cout << "nr = " << real(n) << endl;
    cout << "ni = " << imag(n) << endl;
    cout << "l0 = " << l0 << endl;
    cout << "Q = " << Q << endl;
    cout << "N = " << N << endl;
    cout << "L = " << L << endl;
    cout << "P = " << P << endl;
    cout << "cut_axis = " << cut_axis << endl;
    cout << "cut_value = " << cut_value << endl;
    cout << "xmin = " << xmin << endl;
    cout << "xmax = " << xmax << endl;
    cout << "xpoints = " << xpoints << endl;
    cout << "ymin = " << ymin << endl;
    cout << "ymax = " << ymax << endl;
    cout << "ypoints = " << ypoints << endl;
    cout << "zmin = " << zmin << endl;
    cout << "zmax = " << zmax << endl;
    cout << "zpoints = " << zpoints << endl;
    cout << "file_in = " << file_in << endl;
    cout << "file_out = " << file_out << endl;
    puts("");
    puts("PROGRESS");
    // ===============================================================




    // === F(x,z) FUNCTION ===========================================
    fstream file_mtx( file_in );
    if ( !file_mtx.is_open() )
    {
        cout << "Impossible to open file '" << file_in << "'." << endl;
        puts("The program finished.");
        exit(1);
    }

    // Reading the number of columns and rows of mtx file
    // which also correspond to longitudinal and vertical pixels
    int pxL;
    int pxH;
    file_mtx >> pxH >> pxL;
    
    // If it is not integer, jump to next line.
    while ( !file_mtx.good() )
    {
        file_mtx.clear();
        file_mtx.ignore(INT_MAX, '\n');
        file_mtx >> pxH >> pxL;
    }
    
    int pxHL = pxH*pxL;
    double * F = new double [ pxHL ];
    for ( int i=0; i<pxHL; i++ )
    {
        file_mtx >> F[i];
    }

    file_mtx.close();
    // ===============================================================




    // ==== OUTPUT FILE ==============================================
    file_mtx.open(file_out, fstream::out);
    if ( !file_mtx.is_open() )
    {
        cout << "Impossible to create file '" << file_out << "'." << endl;
        puts("The program finished.");
        exit(1);
    }
    // ===============================================================


    puts("Calculating...");


    // === WAVE VECTORS ==============================================
    int N2M1 = N*2+1;
    complex<double> kk = (M_2PI*n/l0)*(M_2PI*n/l0);
    complex<double> * b = new complex<double> [ N2M1 ]; // beta
    complex<double> * h = new complex<double> [ N2M1 ]; // krho
    
    // In what follows, iq=q+N, q=iq-N, [iq]->q 
    for ( int iq=0; iq<N2M1; iq++ )
    {
        b[iq] = complex<double> ( Q + M_2PI*(iq-N)/L, ( Q + M_2PI*(iq-N)/L )*imag(n)/real(n) );
        h[iq] = sqrt( kk - b[iq]*b[iq] );
    }
    // ===============================================================




    // === COEFFICIENTS Aqps =========================================
    complex<double> * Aqp;
    if ( P < pxH )
    {
        complex<double> n_temp = ( is_nonresistant == 1 ) ? complex<double>( real(n), 0 ) : n; 
        Aqp = new complex<double> [ P*N2M1 ];

        int pxLm1 = pxL-1;

        for ( int ip=0; ip<P; ip++ )
        {
            // Variable that relates pixels in x of input image with number P of LFW
            int pxHip_P = int ( floor(pxH*ip/P) );

            for ( int iq=0; iq<N2M1; iq++ )
            {
                // The following is iteration in j that represents pixels in z. It is an
                // approximation of the Aqp integral by trapezoidal method for equally spaced z.
                // Each j is associated to a point in z (or a pixel in zj of input image).
                complex<double> temp = 0.5*( F[ pxHip_P + pxH*pxLm1 ]*exp( L*Q*imag(n_temp)/real(n_temp) ) + F[ pxHip_P + 0 ] );
                for ( int j=1; j<pxLm1; j++ )
                {
                    temp += F[ pxHip_P + pxH*j ]*exp( complex<double>(0., -M_2PI*(iq-N)*j/pxLm1) )
                            * exp( (Q*imag(n_temp)/real(n_temp))*L*j/pxLm1 );                    
                }
                Aqp[ iq + N2M1*ip ] = temp/( double (pxLm1) );
            }
        }
    }
    else // It is an optmization for pxH < P, Aqp -> A[q][pxh]
    {
        complex<double> n_temp = ( is_nonresistant == 1 ) ? complex<double>( real(n), 0 ) : n; 
        Aqp = new complex<double> [ pxH*N2M1 ];

        int pxLm1 = pxL-1;

        for ( int ipxh=0; ipxh<pxH; ipxh++ )
        {

            for ( int iq=0; iq<N2M1; iq++ )
            {
                // The following is iteration in j that represents pixels in z. It is an
                // approximation of the Aqp integral by trapezoidal method for equally spaced z.
                // Each j is associated to a point in z (or a pixel in zj of input image).
                complex<double> temp = 0.5*( F[ ipxh + pxH*pxLm1 ]*exp( L*Q*imag(n_temp)/real(n_temp) ) + F[ ipxh + 0 ] );
                for ( int j=1; j<pxLm1; j++ )
                {
                    temp += F[ ipxh + pxH*j ]*exp( complex<double>(0., -M_2PI*(iq-N)*j/pxLm1) )
                            * exp( (Q*imag(n_temp)/real(n_temp))*L*j/pxLm1 );                    
                }
                Aqp[ iq + N2M1*ipxh ] = temp/( double (pxLm1) );
            }
        }
    }
    // ===============================================================




    // === DISPLACEMENT IN X OF EACH LFW RELATIVE TO THE ORIGIN ======
    double * x0 = new double [ P ];
    double * x0x0 = new double [ P ];

    double H = L * double(pxH)/double(pxL); 

    // In what follows, ip=p-1, p=ip+1, [ip]->p 
    for ( int ip=0; ip<P; ip++ )
    {
        // Linear and homogeneous distribution in x according to P and H.
        x0[ip] = ( double(ip) )*H/( double(P) );
        x0x0[ip] = x0[ip]*x0[ip];
    }
    // ===============================================================




    // === INTENSITY OF SFW ==========================================
    double * sfwi;
    if ( cut_axis[0] == 'y' || cut_axis[0] == 'Y' )
    {     
        sfwi = new double [ xpoints*zpoints ];
        double dx = (xmax-xmin)/( double(xpoints) );
        double dz = (zmax-zmin)/( double(zpoints) );
        double var_yy = cut_value*cut_value;

        // Longitudinal iteration
        for ( int iz=0; iz<zpoints; iz++ )
        {
            complex<double> var_iz = complex<double> (0., zmin + ( double(iz) ) * dz);

            // Radial iteration
            for ( int ix=0; ix<xpoints; ix++ )
            {
                double var_x = xmin + ( double(ix) ) * dx;
                double var_rhorho = var_x*var_x + var_yy;
                double var_2rho_cosphi = 2.0*var_x;
                complex<double> temp = complex<double> (0.0, 0.0);

                for ( int ip=0; ip<P; ip++ )
                {
                    if ( P < pxH ) // pxH >= P, Aqp -> A[q][p]
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
                        int pxHip_P = int ( floor(pxH*ip/P) );

                        for ( int iq=0; iq<N2M1; iq++ )
                        {
                            temp += Aqp[ iq + N2M1*pxHip_P ]
                                    * bessel::cyl_j0( h[iq]*sqrt(var_rhorho+x0x0[ip]-var_2rho_cosphi*x0[ip]) )
                                    * exp( b[iq]*var_iz );
                        }
                    }
                }
                sfwi[ix + xpoints*iz] = abs(temp)*abs(temp);
            }
        }

    }
    else if ( cut_axis[0] == 'z' || cut_axis[0] == 'Z' )
    {
        sfwi = new double [ xpoints*ypoints ];
        double dx = (xmax-xmin)/( double(xpoints) );
        double dy = (ymax-ymin)/( double(ypoints) );
        complex<double> var_iz = complex<double> (0., cut_value);

        // Transversal iteration
        for ( int iy=0; iy<ypoints; iy++ )
        {
            double var_yy = (ymin + ( double(iy) ) * dy)
                          * (ymin + ( double(iy) ) * dy);

            // Radial iteration 
            for ( int ix=0; ix<xpoints; ix++ )
            {
                double var_x = xmin + ( double(ix) ) * dx;
                double var_rhorho = var_x*var_x + var_yy;
                double var_2rho_cosphi = 2.0*var_x;
                complex<double> temp = complex<double> (0.0, 0.0);

                for ( int ip=0; ip<P; ip++ )
                {
                    if ( P < pxH ) // pxH >= P, Aqp -> A[q][p]
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
                        int pxHip_P = int ( floor(pxH*ip/P) );
                       
                        for ( int iq=0; iq<N2M1; iq++ )
                        {
                            temp += Aqp[ iq + N2M1*pxHip_P ]
                                    * bessel::cyl_j0( h[iq]*sqrt(var_rhorho+x0x0[ip]-var_2rho_cosphi*x0[ip]) )
                                    * exp( b[iq]*var_iz );
                        }
                    }
                }
                sfwi[ix + xpoints*iy] = abs(temp)*abs(temp);
            }
        }
    }
    else if ( cut_axis[0] == 'x' || cut_axis[0] == 'X' )
    {
        sfwi = new double [ ypoints*zpoints ];
        double dy = (ymax-ymin)/( double(ypoints) );
        double dz = (zmax-zmin)/( double(zpoints) );
        double var_xx = cut_value*cut_value;
        double var_2rho_cosphi = 2.0*cut_value;

        // Longitudinal iteration
        for ( int iz=0; iz<zpoints; iz++ )
        {
            complex<double> var_iz = complex<double> (0., zmin + ( double(iz) ) * dz);

            // Transversal iteration
            for ( int iy=0; iy<ypoints; iy++ )
            {
                double var_rhorho = var_xx
                 + (ymin + ( double(iy) ) * dy)*(ymin + ( double(iy) ) * dy); // rho²=x²+y²
                complex<double> temp = complex<double> (0.0, 0.0);

                for ( int ip=0; ip<P; ip++ )
                {
                    if ( P < pxH ) // pxH >= P, Aqp -> A[q][p]
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
                        int pxHip_P = int ( floor(pxH*ip/P) );
                        
                        for ( int iq=0; iq<N2M1; iq++ )
                        {
                            temp += Aqp[ iq + N2M1*pxHip_P ]
                                    * bessel::cyl_j0( h[iq]*sqrt(var_rhorho+x0x0[ip]-var_2rho_cosphi*x0[ip]) )
                                    * exp( b[iq]*var_iz );
                        }
                    }
                }
                sfwi[iy + ypoints*iz] = abs(temp)*abs(temp);
            }
        }
    }
    // ===============================================================




    // === EXPORT INTENSITY ==========================================
    file_mtx << "%%MatrixMarket matrix array real general\n";
    file_mtx << "%nr ni l0 Q N L P H cut_axis cut_value xmin xmax xpoints ";
    file_mtx << "ymin ymax ypoints zmin zmax zpoints is_nonresistant" << endl;  
    file_mtx << "%" << real(n) << " " << imag(n) << " " << l0 << " " << Q;
    file_mtx << " " << N << " " << L << " " << P << " " << H;
    file_mtx << " " << cut_axis << " " << cut_value;
    file_mtx << " " << xmin << " " << xmax << " " << xpoints;
    file_mtx << " " << ymin << " " << ymax << " " << ypoints;
    file_mtx << " " << zmin << " " << zmax << " " << zpoints;
    file_mtx << " " << is_nonresistant << endl;;

    if ( cut_axis[0] == 'y' || cut_axis[0] == 'Y' )
    {        

        file_mtx << xpoints << " " << zpoints << endl;

        for ( int iz=0; iz<zpoints; iz++ )
        {
            for ( int ix=0; ix<xpoints; ix++ )
            {
                file_mtx << setprecision(16) << scientific << "   " << sfwi[ix + xpoints*iz] << endl;
            }
        }
    }
    else if ( cut_axis[0] == 'z' || cut_axis[0] == 'Z' )
    {

        file_mtx << xpoints << " " << ypoints << endl;

        for ( int iy=0; iy<ypoints; iy++ )
        {
            for ( int ix=0; ix<xpoints; ix++ )
            {
                file_mtx << setprecision(16) << scientific << "   " << sfwi[ix + xpoints*iy] << endl;
            }
        }
    }
    else if ( cut_axis[0] == 'x' || cut_axis[0] == 'X' )
    {

        file_mtx << ypoints << " " << zpoints << endl;

        for ( int iz=0; iz<zpoints; iz++ )
        {
            for ( int iy=0; iy<ypoints; iy++ )
            {
                file_mtx << setprecision(16) << scientific << "   " << sfwi[iy + ypoints*iz] << endl;
            }
        }
    }

    file_mtx.close();
    // ===============================================================


    delete[] Aqp;
    delete[] x0;
    delete[] x0x0;
    delete[] b;
    delete[] h;
    delete[] sfwi;
    tempo = time(NULL) - tempo;
    cout << "Elapsed time: " << tempo << " seconds." << endl;

}