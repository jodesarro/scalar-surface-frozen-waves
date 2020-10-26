// REFERENCES
// [1]
// [2] L. A. Ambrosio, “Millimeter-structured nondiffracting surface beams,” J. Opt. Soc. Am. B, vol.
//  36, no. 3, p. 638, Mar. 2019.
// [3] J. O. de Sarro and L. A. Ambrosio, “Constructing Millimeter-structured Surface Beams from
//  Nondiffracting Zeroth-order Bessel Beams in Lossless Media,” in 2019 PhotonIcs & Electromagnetics
//  Research Symposium - Spring (PIERS-Spring), 2019, pp. 283–288.
// [4] M. Zamboni-Rached, “Diffraction-Attenuation resistant beams in absorbing media,” Opt. Express,
//  vol. 14, no. 5, pp. 1804–1809, Mar. 2006.
// [5] M. Zamboni-Rached and M. Mojahedi, “Shaping finite-energy diffraction- and attenuation-resistant
//  beams through Bessel-Gauss–beam superposition,” Phys. Rev. A, vol. 92, no. 4, p. 043839, Oct. 2015.

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

#ifndef M_2PI
#define M_2PI   6.2831853071795864769252867
#endif

#ifndef M_SQRT2
#define M_SQRT2 1.4142135623730950488016887
#endif

using namespace std;

struct InputParameters
{
    complex<double> n;
    double l0, Q, L, H;
    int N, P;
    char* slice_axis_name = new char[1];
    double slice_axis_value, xmin, xmax, ymin, ymax, zmin, zmax;
    int xpoints, ypoints, zpoints;
    string method;
    double w0;
    string file_in, file_out;
    long int start_time;
};

void log_information()
{
    puts("\n| SURFACE-FROZEN-WAVES | 26 OUT 2020 | @JODESARRO |\n");
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

void log_input_parameters( InputParameters parameters )
{
    puts("INPUT PARAMETERS");
    cout << "nr = " << real(parameters.n) << endl;
    cout << "ni = " << imag(parameters.n) << endl;
    cout << "l0 = " << parameters.l0 << endl;
    cout << "Q = " << parameters.Q << endl;
    cout << "N = " << parameters.N << endl;
    cout << "L = " << parameters.L << endl;
    cout << "P = " << parameters.P << endl;
    cout << "slice_axis_name = " << parameters.slice_axis_name[0] << endl;
    cout << "slice_axis_value = " << parameters.slice_axis_value << endl;
    cout << "xmin = " << parameters.xmin << endl;
    cout << "xmax = " << parameters.xmax << endl;
    cout << "xpoints = " << parameters.xpoints << endl;
    cout << "ymin = " << parameters.ymin << endl;
    cout << "ymax = " << parameters.ymax << endl;
    cout << "ypoints = " << parameters.ypoints << endl;
    cout << "zmin = " << parameters.zmin << endl;
    cout << "zmax = " << parameters.zmax << endl;
    cout << "zpoints = " << parameters.zpoints << endl;
    cout << "method = " << parameters.method << endl;
    cout << "w0 = " << parameters.w0 << endl;
    cout << "file_in = " << parameters.file_in << endl;
    cout << "file_out = " << parameters.file_out << endl << endl;
}

void log_progress_calculating()
{
    puts("PROGRESS");
    puts("Calculating...");
}

void calculate_wave_vectors(InputParameters parameters, complex<double> * b, complex<double> * h, int IQMAX)
{
    complex<double> k = M_2PI*parameters.n/parameters.l0;
    complex<double> kk = k*k;
    double k0k0 = ( M_2PI / parameters.l0 ) * ( M_2PI / parameters.l0 );

    // As always in the code, q->[iq] : q=iq-N, 0<=iq<IQMAX.

    if ( parameters.method == "finite_energy" )
    {
        // Eqs (16)-(18) of [5].
        for ( int iq=0; iq<IQMAX; iq++ )
        {
            double q = iq-parameters.N;
            b[iq] = complex<double> ( parameters.Q + M_2PI*q/parameters.L,
                                    imag(k) * ( double(2) - (parameters.Q + M_2PI*q/parameters.L)/real(k) ) );
            h[iq] = M_SQRT2
                    * sqrt( double(1) - (double(1)/real(k)) * (parameters.Q + M_2PI*q/parameters.L) )
                    * abs(k);
        }
    }
    if ( parameters.method == "paraxial_apodized" )
    {
        // Eqs (16)-(18) of [5].
        for ( int iq=0; iq<IQMAX; iq++ )
        {
            double q = iq-parameters.N;
            b[iq] = complex<double> ( parameters.Q + M_2PI*q/parameters.L,
                                    imag(k) * ( double(2) - (parameters.Q + M_2PI*q/parameters.L)/real(k) ) );
            h[iq] = M_SQRT2
                    * sqrt( double(1) - (double(1)/real(k)) * (parameters.Q + M_2PI*q/parameters.L) )
                    * abs(k);
        }
    }
    else if ( parameters.method == "modified" )
    {
        // Eqs (4)-(5) of [5].
        for ( int iq=0; iq<IQMAX; iq++ )
        {
            double q = iq-parameters.N;
            b[iq] = complex<double> ( parameters.Q + M_2PI*q/parameters.L,
                                    k0k0*real(parameters.n)*imag(parameters.n) / (parameters.Q + M_2PI*q/parameters.L) );
            h[iq] = sqrt( (real(parameters.n)*real(parameters.n) - imag(parameters.n)*imag(parameters.n))*k0k0
                          - (parameters.Q + M_2PI*q/parameters.L)*(parameters.Q + M_2PI*q/parameters.L)
                          + (k0k0*real(parameters.n)*imag(parameters.n)/(parameters.Q + M_2PI*q/parameters.L))
                            *(k0k0*real(parameters.n)*imag(parameters.n)/(parameters.Q + M_2PI*q/parameters.L)) );
        }
    }
    else if ( parameters.method == "traditional" || parameters.method == "traditional_nonresistant" )
    {
        // Eqs (2)-(4) of [4].
        for ( int iq=0; iq<IQMAX; iq++ )
        {
            double q = iq-parameters.N;
            b[iq] = complex<double> ( parameters.Q + M_2PI*q/parameters.L,
                                    ( parameters.Q + M_2PI*q/parameters.L ) * imag(parameters.n)/real(parameters.n) );
            h[iq] = sqrt( kk - b[iq]*b[iq] );
        }
    }
}

void calculate_coefficients_Aqp(InputParameters parameters, complex<double> * Aqp, int IQMAX, int IPMAX,
                                double * F, int IIMAX, int IKMAX)
{
    // The following is an approximation of the Aqp integral by trapezoidal method
    //  for equally spaced z. It contains an iteration in ik that represents pixels in z.
    //  Each ik is associated to a point in z.
    // Also, ii is a variable that relates the transversal pixels of the input image with
    //  a number P of LFW once the transversal index is [ii] for F and [ip] for Aqp.

    if ( parameters.method == "finite_energy" )
    {
        // Chapter IV of [5], eq (29) with eqs (16)-(18) and (30).
        complex<double> k = M_2PI*parameters.n/parameters.l0;
        complex<double> h0  = M_SQRT2 * sqrt( double(1) - parameters.Q/real(k) ) * abs(k);
        complex<double> G0 = double (1);
        complex<double> GL = exp( - imag(k)*parameters.L )
                            * exp( - h0*h0 * parameters.L*parameters.L / ( parameters.w0*parameters.w0 * real(k)*real(k) ) );

        for ( int ip=0; ip<IPMAX; ip++ )
        {
            int ii = int ( floor( ip * (IIMAX-1)/(IPMAX-1) ) );

            for ( int iq=0; iq<IQMAX; iq++ )
            {
                double q = iq-parameters.N;

                complex<double> temp = 0.5*( F[ ii + IIMAX*(IKMAX-1) ]/GL + F[ ii + 0 ]/G0 );
                for ( int ik=1; ik<IKMAX-1; ik++ )
                {
                    double z = ik*parameters.L/(IKMAX-1);
                    
                    complex<double> Gz = exp( - imag(k)*z )
                                        * exp( - h0*h0 * z*z / ( parameters.w0*parameters.w0 * real(k)*real(k) ) );

                    temp += ( F[ ii + IIMAX*ik ]/Gz )
                            * exp( complex<double>(0., -M_2PI*q*z/parameters.L) );
                }
                Aqp[ iq + IQMAX*ip ] = temp/( double (IKMAX-1) );
            }
        }
    }
    else if ( parameters.method == "paraxial_apodized" )
    {
        // Chapter III of [5], eq (6) with eqs (16)-(18).
        complex<double> k = M_2PI*parameters.n/parameters.l0;
        double bI0 = imag(k)*( double (2) - parameters.Q/real(k) );

        for ( int ip=0; ip<IPMAX; ip++ )
        {
            int ii = int ( floor( ip * (IIMAX-1)/(IPMAX-1) ) );

            for ( int iq=0; iq<IQMAX; iq++ )
            {
                double q = iq-parameters.N;

                complex<double> temp = 0.5*( ( F[ ii + IIMAX*(IKMAX-1) ] )*exp( bI0*parameters.L ) + F[ ii + 0 ] );
                for ( int ik=1; ik<IKMAX-1; ik++ )
                {
                    double z = ik*parameters.L/(IKMAX-1);
                    
                    temp += ( F[ ii + IIMAX*ik ] )
                            * exp( complex<double>(0., -M_2PI*q*z/parameters.L) )
                            * exp( bI0*z );
                }
                Aqp[ iq + IQMAX*ip ] = temp/( double (IKMAX-1) );
            }
        }
    }
    else if ( parameters.method == "modified" )
    {
        // Chapter II of [5], eq (6) with eqs (4) and (5).
        double k0k0 = ( M_2PI / parameters.l0 ) * ( M_2PI / parameters.l0 );
        double bI0 = k0k0*real(parameters.n)*imag(parameters.n)/parameters.Q;

        for ( int ip=0; ip<IPMAX; ip++ )
        {
            int ii = int ( floor( ip * (IIMAX-1)/(IPMAX-1) ) );

            for ( int iq=0; iq<IQMAX; iq++ )
            {
                double q = iq-parameters.N;

                complex<double> temp = 0.5*( ( F[ ii + IIMAX*(IKMAX-1) ] )*exp( bI0*parameters.L ) + F[ ii + 0 ] );
                for ( int ik=1; ik<IKMAX-1; ik++ )
                {
                    double z = ik*parameters.L/(IKMAX-1);
                    
                    temp += ( F[ ii + IIMAX*ik ] )
                            * exp( complex<double>(0., -M_2PI*q*z/parameters.L) )
                            * exp( bI0*z );
                }
                Aqp[ iq + IQMAX*ip ] = temp/( double (IKMAX-1) );
            }
        }
    }
    else if ( parameters.method == "traditional_nonresistant" )
    {
        // See [3] for what I call a nonresistant SFW.
        // Eq (7) of [4] with eqs (2)-(4).
        for ( int ip=0; ip<IPMAX; ip++ )
        {
            int ii = int ( floor( ip * (IIMAX-1)/(IPMAX-1) ) );

            for ( int iq=0; iq<IQMAX; iq++ )
            {
                double q = iq-parameters.N;

                complex<double> temp = 0.5*( ( F[ ii + IIMAX*(IKMAX-1) ] ) + F[ ii + 0 ] );
                for ( int ik=1; ik<IKMAX-1; ik++ )
                {
                    double z = ik*parameters.L/(IKMAX-1);
                    
                    temp += ( F[ ii + IIMAX*ik ] )
                            * exp( complex<double>(0., -M_2PI*q*z/parameters.L) );
                }
                Aqp[ iq + IQMAX*ip ] = temp/( double (IKMAX-1) );
            }
        }
    }
    else if ( parameters.method == "traditional" )
    {
        // Eq (7) of [4] with eqs (2)-(4).
        double bI0 = parameters.Q * imag(parameters.n)/real(parameters.n);

        for ( int ip=0; ip<IPMAX; ip++ )
        {
            int ii = int ( floor( ip * (IIMAX-1)/(IPMAX-1) ) );

            for ( int iq=0; iq<IQMAX; iq++ )
            {
                double q = iq-parameters.N;

                complex<double> temp = 0.5*( ( F[ ii + IIMAX*(IKMAX-1) ] )*exp( bI0*parameters.L ) + F[ ii + 0 ] );
                for ( int ik=1; ik<IKMAX-1; ik++ )
                {
                    double z = ik*parameters.L/(IKMAX-1);
                    
                    temp += ( F[ ii + IIMAX*ik ] )
                            * exp( complex<double>(0., -M_2PI*q*z/parameters.L) )
                            * exp( bI0*z );
                }
                Aqp[ iq + IQMAX*ip ] = temp/( double (IKMAX-1) );
            }
        }
    }
}

void calculate_displacement_x0p(InputParameters parameters, double * x0, double * x0x0, int IPMAX)
{
    for ( int ip=0; ip<IPMAX; ip++ )
    {
        // Linear and homogeneous distribution in x according to P and H
        x0[ip] = ( double(ip) )*parameters.H/( double(IPMAX) );
        x0x0[ip] = x0[ip]*x0[ip];
    }
}

void calculate_sfw_intensity_xz(InputParameters parameters, double * x0, double * x0x0,
                                complex<double> * b, complex<double> * h,
                                complex<double> * Aqp, int IQMAX, int IPMAX, double * sfwi)
{
    double dx = (parameters.xmax-parameters.xmin)/( double(parameters.xpoints) );
    double dz = (parameters.zmax-parameters.zmin)/( double(parameters.zpoints) );
    double yy = parameters.slice_axis_value*parameters.slice_axis_value;

    if ( parameters.method == "finite_energy" )
    {
        // Eq (20) of [5].
        complex<double> k = M_2PI*parameters.n/parameters.l0;
        for ( int iz=0; iz<parameters.zpoints; iz++ )
        {
            double z = parameters.zmin + ( double(iz) ) * dz;
            complex<double> z_times_i = complex<double> (0., z);
            complex<double> mu = double(1) + complex<double>(0.,1.)*2.0*z/( parameters.w0*parameters.w0 * k );

            for ( int ix=0; ix<parameters.xpoints; ix++ )
            {
                double x = parameters.xmin + ( double(ix) ) * dx;
                double rhorho = x*x + yy;
                double rho_2_cosphi = 2.0*x;
                complex<double> temp = complex<double> (0.0, 0.0);

                for ( int ip=0; ip<IPMAX; ip++ )
                {
                    double w = sqrt( rhorho + x0x0[ip] - rho_2_cosphi*x0[ip] );
                    for ( int iq=0; iq<IQMAX; iq++ )
                    {
                        temp += Aqp[ iq + IQMAX*ip ]
                                * bessel::cyl_j0( h[iq]*w/mu )
                                * exp( b[iq]*z_times_i/mu );
                    }
                    temp *= ( exp(-w*w/(parameters.w0*parameters.w0*mu)) / mu )
                        * exp(-z_times_i*k/mu);
                }
                temp *= exp(z_times_i*k);
                sfwi[ix + parameters.xpoints*iz] = abs(temp)*abs(temp);
            }
        }
    }
    else if ( parameters.method == "paraxial_apodized" )
    {
        // Eq (20) of [5].
        complex<double> k = M_2PI*parameters.n/parameters.l0;
        for ( int iz=0; iz<parameters.zpoints; iz++ )
        {
            double z = parameters.zmin + ( double(iz) ) * dz;
            complex<double> z_times_i = complex<double> (0., z);
            complex<double> mu = double(1) + complex<double>(0.,1.)*2.0*z/( parameters.w0*parameters.w0 * k );

            for ( int ix=0; ix<parameters.xpoints; ix++ )
            {
                double x = parameters.xmin + ( double(ix) ) * dx;
                double rhorho = x*x + yy;
                double rho_2_cosphi = 2.0*x;
                complex<double> temp = complex<double> (0.0, 0.0);

                for ( int ip=0; ip<IPMAX; ip++ )
                {
                    double w = sqrt( rhorho + x0x0[ip] - rho_2_cosphi*x0[ip] );
                    for ( int iq=0; iq<IQMAX; iq++ )
                    {
                        temp += Aqp[ iq + IQMAX*ip ]
                                * bessel::cyl_j0( h[iq]*w/mu )
                                * exp( b[iq]*z_times_i/mu );
                    }
                    temp *= ( exp(-w*w/(parameters.w0*parameters.w0*mu)) / mu )
                        * exp(-z_times_i*k/mu);
                }
                temp *= exp(z_times_i*k);
                sfwi[ix + parameters.xpoints*iz] = abs(temp)*abs(temp);
            }
        }
    }
    else if ( parameters.method == "modified" )
    {
        // Eq (1) of [5].
        for ( int iz=0; iz<parameters.zpoints; iz++ )
        {
            complex<double> z_times_i = complex<double> (0., parameters.zmin + ( double(iz) ) * dz);

            for ( int ix=0; ix<parameters.xpoints; ix++ )
            {
                double x = parameters.xmin + ( double(ix) ) * dx;
                double rhorho = x*x + yy;
                double rho_2_cosphi = 2.0*x;
                complex<double> temp = complex<double> (0.0, 0.0);

                for ( int ip=0; ip<IPMAX; ip++ )
                {
                    double w = sqrt( rhorho + x0x0[ip] - rho_2_cosphi*x0[ip] );
                    for ( int iq=0; iq<IQMAX; iq++ )
                    {
                        temp += Aqp[ iq + IQMAX*ip ]
                                * bessel::cyl_j0( h[iq]*w )
                                * exp( b[iq]*z_times_i );
                    }
                }
                sfwi[ix + parameters.xpoints*iz] = abs(temp)*abs(temp);
            }
        }
    }
    else if ( parameters.method == "traditional" || parameters.method == "traditional_nonresistant" )
    {
        // Eq (6) of [4].
        for ( int iz=0; iz<parameters.zpoints; iz++ )
        {
            complex<double> z_times_i = complex<double> (0., parameters.zmin + ( double(iz) ) * dz);

            for ( int ix=0; ix<parameters.xpoints; ix++ )
            {
                double x = parameters.xmin + ( double(ix) ) * dx;
                double rhorho = x*x + yy;
                double rho_2_cosphi = 2.0*x;
                complex<double> temp = complex<double> (0.0, 0.0);

                for ( int ip=0; ip<IPMAX; ip++ )
                {
                    double w = sqrt( rhorho + x0x0[ip] - rho_2_cosphi*x0[ip] );
                    for ( int iq=0; iq<IQMAX; iq++ )
                    {
                        temp += Aqp[ iq + IQMAX*ip ]
                                * bessel::cyl_j0( h[iq]*w )
                                * exp( b[iq]*z_times_i );
                    }
                }
                sfwi[ix + parameters.xpoints*iz] = abs(temp)*abs(temp);
            }
        }
    }
}

void calculate_sfw_intensity_xy(InputParameters parameters, double * x0, double * x0x0,
                                complex<double> * b, complex<double> * h,
                                complex<double> * Aqp, int IQMAX, int IPMAX, double * sfwi)
{
    double dx = (parameters.xmax-parameters.xmin)/( double(parameters.xpoints) );
    double dy = (parameters.ymax-parameters.ymin)/( double(parameters.ypoints) );
    double z = parameters.slice_axis_value;
    complex<double> z_times_i = complex<double> (0., z);

    if ( parameters.method == "finite_energy" )
    {
        // Eq (20) of [5].
        complex<double> k = M_2PI*parameters.n/parameters.l0;
        complex<double> mu = double(1) + complex<double>(0.,1.)*2.0*z/( parameters.w0*parameters.w0 * k );
        for ( int iy=0; iy<parameters.ypoints; iy++ )
        {
            double yy = (parameters.ymin + ( double(iy) ) * dy)
                            * (parameters.ymin + ( double(iy) ) * dy);

            for ( int ix=0; ix<parameters.xpoints; ix++ )
            {
                double x = parameters.xmin + ( double(ix) ) * dx;
                double rhorho = x*x + yy;
                double rho_2_cosphi = 2.0*x;
                complex<double> temp = complex<double> (0.0, 0.0);

                for ( int ip=0; ip<IPMAX; ip++ )
                {
                    double w = sqrt( rhorho + x0x0[ip] - rho_2_cosphi*x0[ip] );
                    for ( int iq=0; iq<IQMAX; iq++ )
                    {
                        temp += Aqp[ iq + IQMAX*ip ]
                                * bessel::cyl_j0( h[iq]*w/mu )
                                * exp( b[iq]*z_times_i/mu );
                    }
                    temp *= ( exp(-w*w/(parameters.w0*parameters.w0*mu)) / mu )
                        * exp(-z_times_i*k/mu);
                }
                temp *= exp(z_times_i*k);
                sfwi[ix + parameters.xpoints*iy] = abs(temp)*abs(temp);
            }
        }
    }
    else if ( parameters.method == "paraxial_apodized" )
    {
        // Eq (20) of [5].
        complex<double> k = M_2PI*parameters.n/parameters.l0;
        complex<double> mu = double(1) + complex<double>(0.,1.)*2.0*z/( parameters.w0*parameters.w0 * k );
        for ( int iy=0; iy<parameters.ypoints; iy++ )
        {
            double yy = (parameters.ymin + ( double(iy) ) * dy)
                            * (parameters.ymin + ( double(iy) ) * dy);

            for ( int ix=0; ix<parameters.xpoints; ix++ )
            {
                double x = parameters.xmin + ( double(ix) ) * dx;
                double rhorho = x*x + yy;
                double rho_2_cosphi = 2.0*x;
                complex<double> temp = complex<double> (0.0, 0.0);

                for ( int ip=0; ip<IPMAX; ip++ )
                {
                    double w = sqrt( rhorho + x0x0[ip] - rho_2_cosphi*x0[ip] );
                    for ( int iq=0; iq<IQMAX; iq++ )
                    {
                        temp += Aqp[ iq + IQMAX*ip ]
                                * bessel::cyl_j0( h[iq]*w/mu )
                                * exp( b[iq]*z_times_i/mu );
                    }
                    temp *= ( exp(-w*w/(parameters.w0*parameters.w0*mu)) / mu )
                        * exp(-z_times_i*k/mu);
                }
                temp *= exp(z_times_i*k);
                sfwi[ix + parameters.xpoints*iy] = abs(temp)*abs(temp);
            }
        }
    }
    else if ( parameters.method == "modified" )
    {
        // Eq (1) of [5].
        for ( int iy=0; iy<parameters.ypoints; iy++ )
        {
            double yy = (parameters.ymin + ( double(iy) ) * dy)
                            * (parameters.ymin + ( double(iy) ) * dy);

            for ( int ix=0; ix<parameters.xpoints; ix++ )
            {
                double x = parameters.xmin + ( double(ix) ) * dx;
                double rhorho = x*x + yy;
                double rho_2_cosphi = 2.0*x;
                complex<double> temp = complex<double> (0.0, 0.0);

                for ( int ip=0; ip<IPMAX; ip++ )
                {
                    double w = sqrt( rhorho + x0x0[ip] - rho_2_cosphi*x0[ip] );
                    for ( int iq=0; iq<IQMAX; iq++ )
                    {
                        temp += Aqp[ iq + IQMAX*ip ]
                                * bessel::cyl_j0( h[iq]*w )
                                * exp( b[iq]*z_times_i );
                    }
                }
                sfwi[ix + parameters.xpoints*iy] = abs(temp)*abs(temp);
            }
        }
    }
    else if ( parameters.method == "traditional" || parameters.method == "traditional_nonresistant" )
    {
        // Eq (6) of [4].
        for ( int iy=0; iy<parameters.ypoints; iy++ )
        {
            double yy = (parameters.ymin + ( double(iy) ) * dy)
                            * (parameters.ymin + ( double(iy) ) * dy);

            for ( int ix=0; ix<parameters.xpoints; ix++ )
            {
                double x = parameters.xmin + ( double(ix) ) * dx;
                double rhorho = x*x + yy;
                double rho_2_cosphi = 2.0*x;
                complex<double> temp = complex<double> (0.0, 0.0);

                for ( int ip=0; ip<IPMAX; ip++ )
                {
                    double w = sqrt( rhorho + x0x0[ip] - rho_2_cosphi*x0[ip] );
                    for ( int iq=0; iq<IQMAX; iq++ )
                    {
                        temp += Aqp[ iq + IQMAX*ip ]
                                * bessel::cyl_j0( h[iq]*w )
                                * exp( b[iq]*z_times_i );
                    }
                }
                sfwi[ix + parameters.xpoints*iy] = abs(temp)*abs(temp);
            }
        }
    }
}

void calculate_sfw_intensity_yz(InputParameters parameters, double * x0, double * x0x0,
                                complex<double> * b, complex<double> * h,
                                complex<double> * Aqp, int IQMAX, int IPMAX, double * sfwi)
{
    double dy = (parameters.ymax-parameters.ymin)/( double(parameters.ypoints) );
    double dz = (parameters.zmax-parameters.zmin)/( double(parameters.zpoints) );
    double xx = parameters.slice_axis_value*parameters.slice_axis_value;
    double rho_2_cosphi = 2.0*parameters.slice_axis_value;

    if ( parameters.method == "finite_energy" )
    {
        // Eq (20) of [5].
        complex<double> k = M_2PI*parameters.n/parameters.l0;
        for ( int iz=0; iz<parameters.zpoints; iz++ )
        {
            double z = parameters.zmin + ( double(iz) ) * dz;
            complex<double> z_times_i = complex<double> (0., z);
            complex<double> mu = double(1) + complex<double>(0.,1.)*2.0*z/( parameters.w0*parameters.w0 * k );

            for ( int iy=0; iy<parameters.ypoints; iy++ )
            {
                double rhorho = xx
                    + (parameters.ymin + ( double(iy) ) * dy)*(parameters.ymin + ( double(iy) ) * dy);
                complex<double> temp = complex<double> (0.0, 0.0);

                for ( int ip=0; ip<IPMAX; ip++ )
                {
                    double w = sqrt( rhorho + x0x0[ip] - rho_2_cosphi*x0[ip] );
                    for ( int iq=0; iq<IQMAX; iq++ )
                    {
                        temp += Aqp[ iq + IQMAX*ip ]
                                * bessel::cyl_j0( h[iq]*w/mu )
                                * exp( b[iq]*z_times_i/mu );
                    }
                    temp *= ( exp(-w*w/(parameters.w0*parameters.w0*mu)) / mu )
                        * exp(-z_times_i*k/mu);
                }
                temp *= exp(z_times_i*k);
                sfwi[iy + parameters.ypoints*iz] = abs(temp)*abs(temp);
            }
        }
    }
    else if ( parameters.method == "paraxial_apodized" )
    {
        // Eq (20) of [5].
        complex<double> k = M_2PI*parameters.n/parameters.l0;
        for ( int iz=0; iz<parameters.zpoints; iz++ )
        {
            double z = parameters.zmin + ( double(iz) ) * dz;
            complex<double> z_times_i = complex<double> (0., z);
            complex<double> mu = double(1) + complex<double>(0.,1.)*2.0*z/( parameters.w0*parameters.w0 * k );

            for ( int iy=0; iy<parameters.ypoints; iy++ )
            {
                double rhorho = xx
                    + (parameters.ymin + ( double(iy) ) * dy)*(parameters.ymin + ( double(iy) ) * dy);
                complex<double> temp = complex<double> (0.0, 0.0);

                for ( int ip=0; ip<IPMAX; ip++ )
                {
                    double w = sqrt( rhorho + x0x0[ip] - rho_2_cosphi*x0[ip] );
                    for ( int iq=0; iq<IQMAX; iq++ )
                    {
                        temp += Aqp[ iq + IQMAX*ip ]
                                * bessel::cyl_j0( h[iq]*w/mu )
                                * exp( b[iq]*z_times_i/mu );
                    }
                    temp *= ( exp(-w*w/(parameters.w0*parameters.w0*mu)) / mu )
                        * exp(-z_times_i*k/mu);
                }
                temp *= exp(z_times_i*k);
                sfwi[iy + parameters.ypoints*iz] = abs(temp)*abs(temp);
            }
        }
    }
    else if ( parameters.method == "modified" )
    {
        // Eq (1) of [5].
        for ( int iz=0; iz<parameters.zpoints; iz++ )
        {
            complex<double> z_times_i = complex<double> (0., parameters.zmin + ( double(iz) ) * dz);

            for ( int iy=0; iy<parameters.ypoints; iy++ )
            {
                double rhorho = xx
                    + (parameters.ymin + ( double(iy) ) * dy)*(parameters.ymin + ( double(iy) ) * dy);
                complex<double> temp = complex<double> (0.0, 0.0);

                for ( int ip=0; ip<IPMAX; ip++ )
                {
                    double w = sqrt( rhorho + x0x0[ip] - rho_2_cosphi*x0[ip] );
                    for ( int iq=0; iq<IQMAX; iq++ )
                    {
                        temp += Aqp[ iq + IQMAX*ip ]
                                * bessel::cyl_j0( h[iq]*w )
                                * exp( b[iq]*z_times_i );
                    }
                }
                sfwi[iy + parameters.ypoints*iz] = abs(temp)*abs(temp);
            }
        }
    }
    else if ( parameters.method == "traditional" || parameters.method == "traditional_nonresistant" )
    {
        // Eq (6) of [4].
        for ( int iz=0; iz<parameters.zpoints; iz++ )
        {
            complex<double> z_times_i = complex<double> (0., parameters.zmin + ( double(iz) ) * dz);

            for ( int iy=0; iy<parameters.ypoints; iy++ )
            {
                double rhorho = xx
                    + (parameters.ymin + ( double(iy) ) * dy)*(parameters.ymin + ( double(iy) ) * dy);
                complex<double> temp = complex<double> (0.0, 0.0);

                for ( int ip=0; ip<IPMAX; ip++ )
                {
                    double w = sqrt( rhorho + x0x0[ip] - rho_2_cosphi*x0[ip] );
                    for ( int iq=0; iq<IQMAX; iq++ )
                    {
                        temp += Aqp[ iq + IQMAX*ip ]
                                * bessel::cyl_j0( h[iq]*w )
                                * exp( b[iq]*z_times_i );
                    }
                }
                sfwi[iy + parameters.ypoints*iz] = abs(temp)*abs(temp);
            }
        }
    }
}

int main( int argc, char *argv[] )
{

    log_information();

    InputParameters parameters;
    if ( argc != 23 )
    {
        log_wrong_parameters_amount();
        exit(1);
    }
    else
    {
        parameters.start_time = time(NULL);
        parameters.n = complex<double> ( atof(argv[1]), atof(argv[2]) );
        parameters.l0 = atof(argv[3]);
        parameters.Q  = atof(argv[4]);
        parameters.N  = atoi(argv[5]);
        parameters.L  = atof(argv[6]);
        parameters.P  = atoi(argv[7]);
        parameters.slice_axis_name = argv[8];
        parameters.slice_axis_value = atof(argv[9]);
        parameters.xmin = atof(argv[10]);
        parameters.xmax = atof(argv[11]);
        parameters.xpoints = atoi(argv[12]);
        parameters.ymin = atof(argv[13]);
        parameters.ymax = atof(argv[14]);
        parameters.ypoints = atoi(argv[15]);
        parameters.zmin = atof(argv[16]);
        parameters.zmax = atof(argv[17]);
        parameters.zpoints = atoi(argv[18]);
        parameters.method = argv[19];
        parameters.w0 = atof(argv[20]);
        parameters.file_in  = argv[21];
        parameters.file_out = argv[22];
    }

    log_input_parameters( parameters );
 
    log_progress_calculating();

    // Read MTX file to function F[ii][ik]
    // ii : 0<=ii<IIMAX
    // ik : 0<=ik<IKMAX
    int IIMAX, IKMAX;
    mtx_get_dimensions( parameters.file_in, IIMAX, IKMAX );
    double * F = new double [ IIMAX*IKMAX ];
    mtx_read_to_array( parameters.file_in, F );

    // Resetting gaussian waist of the Bessel-Gauss beam
    if ( parameters.method != "finite_energy" && parameters.method != "paraxial_apodized" )
    {
        parameters.w0 = 0;
    }

    // Values for q and p indexes of FWs
    // q->[iq] : q=iq-N, 0<=iq<IQMAX
    // p->[ip] : p=ip+1, 0<=ip<IPMAX
    int IQMAX = parameters.N*2+1;
    int IPMAX = parameters.P;

    // Wave vectors
    complex<double> * b = new complex<double> [ IQMAX ];
    complex<double> * h = new complex<double> [ IQMAX ];
    calculate_wave_vectors(parameters, b, h, IQMAX);
    
    // Coefficients Aqp
    complex<double> * Aqp = new complex<double> [ IQMAX*IPMAX ];
    calculate_coefficients_Aqp(parameters, Aqp, IQMAX, IPMAX, F, IIMAX, IKMAX);
    delete[] F;

    // Displacement in x of each LFW relative to the origin
    parameters.H = parameters.L * double(IIMAX)/double(IKMAX); 
    double * x0 = new double [ IPMAX ];
    double * x0x0 = new double [ IPMAX ];
    calculate_displacement_x0p(parameters, x0, x0x0, IPMAX);

    // Intensity of SFW
    double * sfwi;
    int mtx_rows, mtx_columns;
    if ( parameters.slice_axis_name[0] == 'y' )
    {     
        parameters.ymin = 0;
        parameters.ymax = 0;
        parameters.ypoints = 0;
        mtx_rows = parameters.xpoints;
        mtx_columns = parameters.zpoints;
        sfwi = new double [ mtx_rows*mtx_columns ];
        calculate_sfw_intensity_xz(parameters, x0, x0x0, b, h, Aqp, IQMAX, IPMAX, sfwi);
    }
    else if ( parameters.slice_axis_name[0] == 'z' )
    {
        parameters.zmin = 0;
        parameters.zmax = 0;
        parameters.zpoints = 0;
        mtx_rows = parameters.xpoints;
        mtx_columns = parameters.ypoints;
        sfwi = new double [ mtx_rows*mtx_columns ];
        calculate_sfw_intensity_xy(parameters, x0, x0x0, b, h, Aqp, IQMAX, IPMAX, sfwi);
    }
    else if ( parameters.slice_axis_name[0] == 'x' )
    {
        parameters.xmin = 0;
        parameters.xmax = 0;
        parameters.xpoints = 0;
        mtx_rows = parameters.ypoints;
        mtx_columns = parameters.zpoints;
        sfwi = new double [ mtx_rows*mtx_columns ];
        calculate_sfw_intensity_yz(parameters, x0, x0x0, b, h, Aqp, IQMAX, IPMAX, sfwi);
    }
    
    delete[] Aqp;
    delete[] x0;
    delete[] x0x0;
    delete[] b;
    delete[] h;

    // Write to .mtx file
    string mtx_comments;
    mtx_comments = "nr ni l0 Q N L P H slice_axis_name slice_axis_value "
                    "xmin xmax xpoints ymin ymax ypoints zmin zmax zpoints "
                    "method w0 elapsed_time\n"
                    "%" + to_string(real(parameters.n)) + " " + to_string(imag(parameters.n)) + " "
                    + to_string(parameters.l0) + " " + to_string(parameters.Q) + " "
                    + to_string(parameters.N) + " " + to_string(parameters.L) + " "
                    + to_string(parameters.P) + " " + to_string(parameters.H) + " "
                    + parameters.slice_axis_name[0] + " " + to_string(parameters.slice_axis_value) + " "
                    + to_string(parameters.xmin) + " " + to_string(parameters.xmax) + " "
                    + to_string(parameters.xpoints) + " " + to_string(parameters.ymin) + " "
                    + to_string(parameters.ymax) + " " + to_string(parameters.ypoints) + " "
                    + to_string(parameters.zmin) + " " + to_string(parameters.zmax) + " "
                    + to_string(parameters.zpoints) + " " + parameters.method + " "
                    + to_string(parameters.w0) + " " + to_string( time(NULL) - parameters.start_time );
    mtx_write_from_array( parameters.file_out, mtx_comments, mtx_rows, mtx_columns, sfwi );

    delete[] sfwi;
    
    cout << "Done.\nElapsed time: " << ( time(NULL) - parameters.start_time ) << " seconds." << endl;
}