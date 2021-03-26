// REFERENCES *&#946;&#x035E;*<sub>*p*</sub>
// [1] J. O. de Sarro, L. A. Ambrosio, "Surface beams resistant to diffraction and attenuation and
//  structured at the millimeter scale," J. Opt. Soc. Am. B 38, 677-684 (2021).
// [3] J. O. de Sarro and L. A. Ambrosio, "Constructing Millimeter-structured Surface Beams from
//  Nondiffracting Zeroth-order Bessel Beams in Lossless Media," in 2019 PhotonIcs & Electromagnetics
//  Research Symposium - Spring (PIERS-Spring), 2019, pp. 283–288.
// [5] M. Zamboni-Rached and M. Mojahedi, "Shaping finite-energy diffraction- and attenuation-resistant
//  beams through Bessel-Gauss–beam superposition," Phys. Rev. A, vol. 92, no. 4, p. 043839, Oct. 2015.

#include <iostream>
#include <complex>
#include <iomanip>
#include <ctime>
#include <regex>
#include <future>
#include <vector>
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

struct Parameters
{
    // General parameters of FWs
    string method, polarization;
    double nr, ni;
    complex<double> n;
    double l0, k0;
    complex <double> k;
    int N, P;
    double dx0, lfw_spot_radius_max;
    double * L;
    double * Q;
    double * w;
    double * x0;

    // Calculation parameters
    double xmin, xmax, ymin, ymax, zmin, zmax;
    int xpoints, ypoints, zpoints;
    bool is_async;
};

void log_information()
{
    puts("| SURFACE-FROZEN-WAVES | 25 MAR 2021 | @JODESARRO |\n");
    puts("DESCRIPTION");
    puts("A C++ routine to evaluate the fields of zeroth order "
         "surface frozen waves. Documentation and latest release available "
         "on 'https://github.com/jodesarro/surface-frozen-waves'.\n");
}

void log_parameters( Parameters pm, int IQMAX, int IPMAX, ostream& stream )
{
    stream << setprecision(15);
    stream << "PARAMETERS" << endl;
    stream << "method = " << pm.method << endl;
    stream << "polarization = " << pm.polarization << endl;
    stream << "n = " << pm.n << endl;
    stream << "l0 = " << pm.l0 << endl;
    stream << "k0 = " << pm.k0 << endl;
    stream << "k = " << pm.k << endl;
    stream << "N = " << pm.N << endl;
    stream << "P = " << pm.P << endl;
    stream << "lfw_spot_radius_max = " << pm.lfw_spot_radius_max << endl;
    stream << "dx0 = " << pm.dx0 << endl;

    stream << endl;

    stream << "L = {" << pm.L[0];
    for ( int ip=1; ip<IPMAX; ip++ )
    {
        stream << "," << pm.L[ip];
    }
    stream << "}\n" << endl;
   
    stream << "Q = {" << pm.Q[0];
    for ( int ip=1; ip<IPMAX; ip++ )
    {
        stream << "," << pm.Q[ip];
    }
    stream << "}\n" << endl;

    stream << "Nmax = {" << floor( pm.L[0] * min( pm.k0*pm.nr - pm.Q[0] , pm.Q[0] ) / M_2PI );
    for ( int ip=1; ip<IPMAX; ip++ )
    {
        stream << "," << floor( pm.L[ip] * min( pm.k0*pm.nr - pm.Q[ip] , pm.Q[ip] ) / M_2PI );
    }
    stream << "}\n" << endl;

    stream << "w = {" << pm.w[0];
    for ( int ip=1; ip<IPMAX; ip++ )
    {
        stream << "," << pm.w[ip];
    }
    stream << "}\n" << endl;

    stream << "delta = {" << ( double(2) ) * M_2PI * ( double(pm.N) ) /( pm.L[0] * pm.Q[0] );
    for ( int ip=1; ip<IPMAX; ip++ )
    {
        stream << "," << ( double(2) ) * M_2PI * ( double(pm.N) ) /( pm.L[ip] * pm.Q[ip] );
    }
    stream << "}\n" << endl;
    
    stream << "DATA RANGE PARAMETERS" << endl;
    stream << "xmin = " << pm.xmin << endl;
    stream << "xmax = " << pm.xmax << endl;
    stream << "xpoints = " << pm.xpoints << endl;
    stream << "ymin = " << pm.ymin << endl;
    stream << "ymax = " << pm.ymax << endl;
    stream << "ypoints = " << pm.ypoints << endl;
    stream << "zmin = " << pm.zmin << endl;
    stream << "zmax = " << pm.zmax << endl;
    stream << "zpoints = " << pm.zpoints << endl;

    stream << endl;
}

void log_progress( )
{
    puts("PROGRESS");
    puts("Calculating the fields...");
}

void calculate_wave_numbers( Parameters pm,
                                int IQMAX, int IPMAX,
                                complex<double> * h, complex<double> * b )
{
    // Just to remember:
    //  q->[iq] : q=iq-N, 0<=iq<IQMAX
    //  p->[ip] : p=ip+1, 0<=ip<IPMAX

    if ( pm.method == "finite_energy" )
    {
        // Eqs (16)-(18) of [5].
        for ( int ip=0; ip<IPMAX; ip++ )
        {
            for ( int iq=0; iq<IQMAX; iq++ )
            {
                double q = double( iq - pm.N );
                double bqp_r = pm.Q[ip] + M_2PI*q/pm.L[ip];
                b[iq + IQMAX*ip] = complex<double> (
                                                    bqp_r,
                                                    imag(pm.k)*( double(2) - bqp_r/real(pm.k) )
                                                   );
                h[iq + IQMAX*ip] = M_SQRT2
                                    * sqrt( double(1) - bqp_r/real(pm.k) )
                                    * abs(pm.k);
            }
        }
    }
    else if ( pm.method == "resistant_realh" )
    {
        // Eqs (4)-(5) of [5].
        double k0k0 = pm.k0 * pm.k0;
        for ( int ip=0; ip<IPMAX; ip++ )
        {
            for ( int iq=0; iq<IQMAX; iq++ )
            {
                double q = double( iq - pm.N );
                double bqp_r = pm.Q[ip] + M_2PI*q/pm.L[ip];
                b[iq + IQMAX*ip] = complex<double> (
                                                    bqp_r,
                                                    k0k0*pm.nr*pm.ni / bqp_r
                                                   );
                h[iq + IQMAX*ip] = sqrt(
                                        ( pm.nr*pm.nr - pm.ni*pm.ni )*k0k0
                                        - ( bqp_r * bqp_r )
                                        + ( k0k0*pm.nr*pm.ni/bqp_r )*( k0k0*pm.nr*pm.ni/bqp_r )
                                       );
            }
        }
    }
    else if ( pm.method == "resistant" || pm.method == "nonresistant" )
    {
        // Eq (3) of [1].
        complex<double> kk = pm.k * pm.k;
        for ( int ip=0; ip<IPMAX; ip++ )
        {
            for ( int iq=0; iq<IQMAX; iq++ )
            {
                double q = double( iq - pm.N );
                double bqp_r = pm.Q[ip] + M_2PI*q/pm.L[ip];
                b[iq + IQMAX*ip] = complex<double> (
                                                    bqp_r,
                                                    bqp_r * pm.ni/pm.nr
                                                   );
                h[iq + IQMAX*ip] = sqrt(
                                        kk - b[iq + IQMAX*ip]*b[iq + IQMAX*ip]
                                       );
            }
        }
    }
}

void calculate_coefficients_A( Parameters pm, complex<double> * h, complex<double> * b,
                                int IQMAX, int IPMAX, double * F, int IIMAX, int IKMAX,
                                complex<double> * A )
{
    // The following is an approximation of the A_qp integral by trapezoidal method
    //  for equally spaced z. It contains an iteration in ik that represents pixels in z.
    //  Each ik is associated to a point in z.
    // Also, ii is a variable that relates the transversal pixels of the SIP with
    //  a number P of LFW once the transversal index is [ii] for F but [ip] for A.

    if ( pm.method == "finite_energy" )
    {
        // Chapter IV of [5], eq (29) with eqs (16)-(18) and (30).

        double kk_r = real(pm.k) * real(pm.k);
        for ( int ip=0; ip<IPMAX; ip++ )
        {
            int ii = int ( floor( ( double(ip) ) * ( double(IIMAX-1) )/( double(IPMAX-1) ) ) );
            double ww = pm.w[ip]*pm.w[ip];
            complex<double> h0h0 = h[pm.N + IQMAX*ip] * h[pm.N + IQMAX*ip];
            complex<double> G0 = double(1);
            complex<double> GL = exp(   - imag(pm.k)*pm.L[ip]
                                        - h0h0*pm.L[ip]*pm.L[ip]/( ww*kk_r )
                                     );

            for ( int iq=0; iq<IQMAX; iq++ )
            {
                double q = double( iq - pm.N );
                complex<double> aux = (
                                            F[ ii + IIMAX*(IKMAX-1) ]/GL
                                            + F[ ii + 0 ]/G0
                                        )/( double(2) );
                for ( int ik=1; ik<IKMAX-1; ik++ )
                {
                    double z = ( double(ik) ) * pm.L[ip] / ( double(IKMAX-1) );
                    complex<double> Gz = exp( - imag(pm.k)*z
                                              - h0h0*z*z/( ww*kk_r )
                                             );

                    aux += ( F[ ii + IIMAX*ik ]/Gz )
                            * exp( complex<double>(0, -M_2PI*q*z/pm.L[ip]) );
                }
                A[ iq + IQMAX*ip ] = aux/( double(IKMAX-1) );
            }
        }
    }
    else if ( pm.method == "nonresistant" || pm.method == "resistant" || pm.method == "resistant_realh" )
    {
        // Eq (4) of [1].
    
        for ( int ip=0; ip<IPMAX; ip++ )
        {
            int ii = int ( floor( ( double(ip) ) * ( double(IIMAX-1) )/( double(IPMAX-1) ) ) );

            // See [3] or section 3.A of [1] for what I call a nonresistant SFW.
            double beta0_i = imag(b[pm.N + IQMAX*ip]);
            if ( pm.method == "nonresistant" )
            {
                beta0_i = 0.;
            }

            for ( int iq=0; iq<IQMAX; iq++ )
            {
                double q = double( iq - pm.N );
                complex<double> aux = (
                                            ( F[ ii + IIMAX*(IKMAX-1) ] )*exp( beta0_i * pm.L[ip] )
                                            + F[ ii + 0 ]
                                        )/( double(2) );
                for ( int ik=1; ik<IKMAX-1; ik++ )
                {
                    double z = ( double(ik) ) * pm.L[ip] / ( double(IKMAX-1) );
                    aux += ( F[ ii + IIMAX*ik ] )
                            * exp( complex<double>(0, -M_2PI*q*z/pm.L[ip]) )
                            * exp( beta0_i * z );
                }
                A[ iq + IQMAX*ip ] = aux/( double(IKMAX-1) );
            }
        }
    }
}

void fields_finite_energy_in_xy( int iz, double dx, double dy, double dz,
                                Parameters pm, complex<double> * A,
                                complex<double> * h, complex<double> * b, int IQMAX, int IPMAX,
                                complex<double> * field_Ex, complex<double> * field_Ey, complex<double> * field_Ez )
{
    // Eq (20) of [5].

    double z = pm.zmin + ( double(iz) ) * dz;
    complex<double> z_times_i = complex<double> (0, z);

    for ( int iy=0; iy<pm.ypoints; iy++ )
    {
        double y = pm.ymin + ( double(iy) ) * dy;
        double yy = y * y;

        for ( int ix=0; ix<pm.xpoints; ix++ )
        {

            double x = pm.xmin + ( double(ix) ) * dx;
            double xx = x * x;
            
            double rhorho = xx + yy;
            double rho_2_cosphi = ( double(2) ) * x;

            complex<double> psi_parallel = complex<double> (0, 0);
            complex<double> psi_perpendicular = complex<double> (0, 0);

            for ( int ip=0; ip<IPMAX; ip++ )
            {
                double ww = pm.w[ip]*pm.w[ip];
                complex<double> psi = complex<double> (0, 0);
                complex<double> mu = double(1) + ( complex<double>(0,2) )*z/( ww * pm.k );
                double rhoprhop = rhorho + pm.x0[ip]*pm.x0[ip] - rho_2_cosphi*pm.x0[ip];
                double rhop = sqrt( rhoprhop );
                for ( int iq=0; iq<IQMAX; iq++ )
                {
                    psi += A[ iq + IQMAX*ip ]
                            * bessel::cyl_j0( h[iq + IQMAX*ip]*rhop/mu )
                            * exp(
                                    + b[iq + IQMAX*ip]*z_times_i/mu
                                    - z_times_i*pm.k/mu
                                    + z_times_i*pm.k
                                    - rhoprhop/(ww*mu) 
                                )
                            / mu;
                }

                // To distinguish parallel and perpendicular polarizations
                if ( ip % 2 ) // ip odd -> p even
                {
                    psi_perpendicular += psi;
                }
                else // ip even -> p odd
                {
                    psi_parallel += psi;
                }
            }

            if ( pm.polarization == "scalar" )
            {
                field_Ex[ix + pm.xpoints*(iy + pm.ypoints*iz)] = psi_parallel + psi_perpendicular;
                field_Ey[ix + pm.xpoints*(iy + pm.ypoints*iz)] = complex<double>(0,0);
                field_Ez[ix + pm.xpoints*(iy + pm.ypoints*iz)] = complex<double>(0,0);
            }
            if ( pm.polarization == "linear" )
            {
                field_Ex[ix + pm.xpoints*(iy + pm.ypoints*iz)] = psi_parallel + psi_perpendicular;
                field_Ey[ix + pm.xpoints*(iy + pm.ypoints*iz)] = complex<double>(0,0);
                field_Ez[ix + pm.xpoints*(iy + pm.ypoints*iz)] = complex<double>(0,0); // Paraxial beam
            }
            else if ( pm.polarization == "linear_crossed" )
            {
                field_Ex[ix + pm.xpoints*(iy + pm.ypoints*iz)] = psi_parallel;
                field_Ey[ix + pm.xpoints*(iy + pm.ypoints*iz)] = psi_perpendicular;
                field_Ez[ix + pm.xpoints*(iy + pm.ypoints*iz)] = complex<double>(0,0); // Paraxial beam
            }
        }
    }
}

void fields_resistant_in_xy( int iz, double dx, double dy, double dz,
                            Parameters pm, complex<double> * A,
                            complex<double> * h, complex<double> * b,
                            int IQMAX, int IPMAX,
                            complex<double> * field_Ex, complex<double> * field_Ey, complex<double> * field_Ez )
{
    // Eq (2) of [1].

    double z = pm.zmin + ( double(iz) ) * dz;
    complex<double> z_times_i = complex<double> (0, z);

    for ( int iy=0; iy<pm.ypoints; iy++ )
    {
        double y = pm.ymin + ( double(iy) ) * dy;
        double yy = y * y;
        
        for ( int ix=0; ix<pm.xpoints; ix++ )
        {
            double x = pm.xmin + ( double(ix) ) * dx;
            double xx = x * x;
            
            double rhorho = xx + yy;
            double rho_2_cosphi = ( double(2) ) * x;

            complex<double> psi_parallel = complex<double> (0, 0);
            complex<double> psi_perpendicular = complex<double> (0, 0);
            
            for ( int ip=0; ip<IPMAX; ip++ )
            {
                complex<double> psi = complex<double> (0, 0);
                double rhop = sqrt( rhorho + pm.x0[ip]*pm.x0[ip] - rho_2_cosphi*pm.x0[ip] );
                for ( int iq=0; iq<IQMAX; iq++ )
                {
                    psi += A[ iq + IQMAX*ip ]
                            * bessel::cyl_j0( h[iq + IQMAX*ip]*rhop )
                            * exp( b[iq + IQMAX*ip]*z_times_i );
                }

                // To distinguish parallel and perpendicular polarizations
                if ( ip % 2 ) // ip odd -> p even
                {
                    psi_perpendicular += psi;
                }
                else // ip even -> p odd
                {
                    psi_parallel += psi;
                }

            }

            if ( pm.polarization == "scalar" )
            {
                field_Ex[ix + pm.xpoints*(iy + pm.ypoints*iz)] = psi_parallel + psi_perpendicular;
                field_Ey[ix + pm.xpoints*(iy + pm.ypoints*iz)] = complex<double>(0,0);
                field_Ez[ix + pm.xpoints*(iy + pm.ypoints*iz)] = complex<double>(0,0);
            }
            else if ( pm.polarization == "linear" )
            {
                field_Ex[ix + pm.xpoints*(iy + pm.ypoints*iz)] = psi_parallel + psi_perpendicular;
                field_Ey[ix + pm.xpoints*(iy + pm.ypoints*iz)] = complex<double>(0,0);
                field_Ez[ix + pm.xpoints*(iy + pm.ypoints*iz)] = complex<double>(0,0); // Paraxial beam
            }
            else if ( pm.polarization == "linear_crossed" )
            {
                field_Ex[ix + pm.xpoints*(iy + pm.ypoints*iz)] = psi_parallel;
                field_Ey[ix + pm.xpoints*(iy + pm.ypoints*iz)] = psi_perpendicular;
                field_Ez[ix + pm.xpoints*(iy + pm.ypoints*iz)] = complex<double>(0,0); // Paraxial beam
            }
        }
    }
}

void calculate_field( Parameters pm, complex<double> * A,
                        complex<double> * h, complex<double> * b, int IQMAX, int IPMAX,
                        complex<double> * field_Ex, complex<double> * field_Ey, complex<double> * field_Ez )
{
    double dx = (pm.xmax-pm.xmin)/( double(pm.xpoints-1) );
    double dy = (pm.ymax-pm.ymin)/( double(pm.ypoints-1) );
    double dz = (pm.zmax-pm.zmin)/( double(pm.zpoints-1) );

    if ( pm.is_async )
    {
        std::vector<std::future<void>> futures;
        if ( pm.method == "finite_energy" )
        {
            for ( int iz=0; iz<pm.zpoints; iz++ )
            {
                futures.push_back ( async( launch::async, fields_finite_energy_in_xy, iz, dx, dy, dz, pm, A, h, b, IQMAX, IPMAX, field_Ex, field_Ey, field_Ez ) );
            }       
        }
        else if ( pm.method == "resistant" || pm.method == "nonresistant" || pm.method == "resistant_realh" )
        {
            for ( int iz=0; iz<pm.zpoints; iz++ )
            {
                futures.push_back ( async( launch::async, fields_resistant_in_xy, iz, dx, dy, dz, pm, A, h, b, IQMAX, IPMAX, field_Ex, field_Ey, field_Ez ) );
            } 
        }
        for(auto &e : futures)
        {
            e.get();
        }
    }
    else
    {
        if ( pm.method == "finite_energy" )
        {
            for ( int iz=0; iz<pm.zpoints; iz++ )
            {
                fields_finite_energy_in_xy( iz, dx, dy, dz, pm, A, h, b, IQMAX, IPMAX, field_Ex, field_Ey, field_Ez );
            }       
        }
        else if ( pm.method == "resistant" || pm.method == "nonresistant" || pm.method == "resistant_realh" )
        {
            for ( int iz=0; iz<pm.zpoints; iz++ )
            {
                fields_resistant_in_xy( iz, dx, dy, dz, pm, A, h, b, IQMAX, IPMAX, field_Ex, field_Ey, field_Ez );
            } 
        }
    }
}

template<typename T>
void m_3d_write_from_complex_array( string m_path, string m_comment, int IMAX, int JMAX, int KMAX, complex<T> * m_array )
{    
    fstream m_file;
    m_file.open( m_path, fstream::out );
    m_file << "(*" << m_comment << "*)" << endl;
    m_file << "{";
    for ( int i=0; i<IMAX; i++ )
    {
        m_file << "{";
        for ( int j=0; j<JMAX; j++ )
        {
            m_file << "{";
            for ( int k=0; k<KMAX; k++ )
            {
                stringstream number;
                number << scientific << setprecision(16) << showpos ;
                number << real(m_array[i + IMAX*(j + JMAX*k)]) << imag(m_array[i + IMAX*(j + JMAX*k)]) << "*I";
                m_file << regex_replace(number.str(), regex{"e"}, "*^");
                ( k == KMAX-1 ) ? ( m_file << "}") : ( m_file << ",") ;
            }
            ( j == JMAX-1 ) ? ( m_file << "}") : ( m_file << ",") ;
        }
        ( i == IMAX-1 ) ? ( m_file << "}") : ( m_file << ",") ;
    }
    m_file.close();
}

int main()
{
    long int start_time = time( nullptr );

    log_information();

    Parameters pm;

    // General parameters of FWs
    pm.method = "resistant";
    pm.polarization = "scalar";
    pm.nr = 1.4;
    pm.ni = 256e-8;
    pm.n = complex<double> ( pm.nr , pm.ni );
    pm.l0 = 632e-9;
    pm.k0 = M_2PI / pm.l0;
    pm.k = pm.k0 * pm.n;
    pm.N  = 15;
    pm.P  = 472;

    // Maximum values for q and p indexes of FWs
    //  q->[iq] : q=iq-N, 0<=iq<IQMAX
    //  p->[ip] : p=ip+1, 0<=ip<IPMAX
    int IQMAX = pm.N*2+1;
    int IPMAX = pm.P;

    // Parameter Lp
    pm.L = new double [ IPMAX ];
    for ( int ip=0; ip<IPMAX; ip++ )
    {
        pm.L[ip] = 0.06;
    }

    // Parameter Qp
    pm.Q = new double [ IPMAX ];
    for ( int ip=0; ip<IPMAX; ip++ )
    {
        pm.Q[ip] = 0.999 * pm.nr * pm.k0;
        //double r0 = 0.000050;
        //pm.Q[ip] = ( double(1) - 0.5*((2.405)/(r0*abs(pm.k)))*((2.405)/(r0*abs(pm.k))) ) * real(pm.k);
    }

    // Wavenumbers
    complex<double> * h = new complex<double> [ IQMAX*IPMAX ];
    complex<double> * b = new complex<double> [ IQMAX*IPMAX ];
    calculate_wave_numbers( pm, IQMAX, IPMAX, h, b );

    // Parameter wp
    pm.w = new double [ IPMAX ];
    for ( int ip=0; ip<IPMAX; ip++ )
    {
        pm.w[ip] = ( pm.L[ip] * real( h[pm.N + IQMAX*ip] ) )/( ( double(2) )*real(pm.k) );
    }

    // Finding the maximum spot among all P LFWs to eventually use in dx0
    pm.lfw_spot_radius_max = double(0);
    for ( int ip=0; ip<IPMAX; ip++ )
    {
        double aux = 2.405 / real( h[pm.N + IQMAX*ip] );
        if ( aux > pm.lfw_spot_radius_max )
        {
            pm.lfw_spot_radius_max = aux;
        }
    }

    // Separation in x
    pm.dx0 = 0.00003454310980103169;
    //pm.dx0 = 2.0 * pm.lfw_spot_radius_max;
    //pm.dx0 = 128e-6;

    // Displacement in x of each LFW relative to the origin
    pm.x0 = new double [ IPMAX ];
    for ( int ip=0; ip<IPMAX; ip++ )
    {
        pm.x0[ip] = ( double(ip) ) * pm.dx0;
    }

    // Data calculation parameters
    pm.xmin = 0.;
    pm.xmax = ( double(pm.P) ) * pm.dx0;
    pm.xpoints = 208;
    pm.ymin = 0.;
    pm.ymax = 0.5*pm.xmax;
    pm.ypoints = 104;
    pm.zmin = 0.;
    pm.zmax = pm.L[0];
    pm.zpoints = 276;

    // Asynchronously calculation
    pm.is_async = false;

    // There is no need to change values in all of the following

    log_parameters( pm, IQMAX, IPMAX, cout );

    ofstream file_parameters("parameters.txt");
    log_parameters( pm, IQMAX, IPMAX, file_parameters );
    file_parameters.close();

    log_progress();

    // Read MTX SIP file to function F[ii][ik]
    // ii : 0<=ii<IIMAX
    // ik : 0<=ik<IKMAX
    int IIMAX, IKMAX;
    mtx_get_dimensions( "sip.mtx", IIMAX, IKMAX );
    double * F = new double [ IIMAX*IKMAX ];
    mtx_read_to_array( "sip.mtx", F );

    // Coefficients A
    complex<double> * A = new complex<double> [ IQMAX*IPMAX ];
    calculate_coefficients_A( pm, h, b, IQMAX, IPMAX, F, IIMAX, IKMAX, A );
    delete[] F;

    // Field calculation
    complex<double> * field_Ex = new complex<double> [ pm.xpoints*pm.ypoints*pm.zpoints ];
    complex<double> * field_Ey = new complex<double> [ pm.xpoints*pm.ypoints*pm.zpoints ];
    complex<double> * field_Ez = new complex<double> [ pm.xpoints*pm.ypoints*pm.zpoints ];
    calculate_field( pm, A, h, b, IQMAX, IPMAX, field_Ex, field_Ey, field_Ez );

    // Write psi to .m file
    m_3d_write_from_complex_array( "Ex.m", " SFW: component x of the electric field ", pm.xpoints, pm.ypoints, pm.zpoints, field_Ex );
    m_3d_write_from_complex_array( "Ey.m", " SFW: component y of the electric field ", pm.xpoints, pm.ypoints, pm.zpoints, field_Ey );
    m_3d_write_from_complex_array( "Ez.m", " SFW: component z of the electric field ", pm.xpoints, pm.ypoints, pm.zpoints, field_Ez );

    cout << "Done.\nElapsed time: " << ( time(nullptr) - start_time ) << " seconds." << endl;

    system("pause");
}
