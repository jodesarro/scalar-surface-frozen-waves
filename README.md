# Scalar surface frozen waves
A C++ routine to evaluate the scalar field of zeroth order surface frozen waves (SFWs) [[1-5](#references)].

The ordinary scalar SFW is composed by linear frozen waves (LFWs) equally spaced, with the same parameter *N* but distincts *Q*<sub>*p*</sub> and *L*<sub>*p*</sub>.

## Usage
After compiling, the program must be executed with the following parameters passed as arguments, '*nr* *ni* *l0* *file_Q* *N* *file_L* *P* *dx0* *slice_axis_name* *slice_axis_value* *xmin* *xmax* *xpoints* *ymin* *ymax* *ypoints* *zmin* *zmax* *zpoints* *method* *w0* *file_bbar* *file_sip* *file_psi*', where:

- ***nr***: real part of the refractive index;
- ***ni***: imaginary part of the refractive index;
- ***l0***: wavelength in vacuum;
- ***file_Q***: file path to load the parameter *Q*<sub>*p*</sub> in <a href="https://math.nist.gov/MatrixMarket/formats.html">Matrix Market Array File Format (MTX)</a>;
- ***N***: number related to Bessel beams, resulting in a total of 2*N*+1 Bessel beams;
- ***file_L***: file path to load the parameter *L*<sub>*p*</sub> in <a href="https://math.nist.gov/MatrixMarket/formats.html">MTX</a>;
- ***P***: number of LFWs;
- ***dx0***: separation in *x* between consecutives LFWs;
- ***slice_axis_name***: the axis (*x*, *y* or *z*) where the 2D-plot occur perpendicular to it;
- ***slice_axis_value***: the point of the *slice_axis_name* where the 2D-plot occur;
- ***xmin***: minimum *x* to be calculated (it will be ignored if *slice_axis_name* = *x*);
- ***xmax***: maximum *x* to be calculated (it will be ignored if *slice_axis_name* = *x*);
- ***xpoints***: number of points in *x* (it will be ignored if *slice_axis_name* = *x*);
- ***ymin***: minimum *y* to be calculated (it will be ignored if *slice_axis_name* = *y*);
- ***ymax***: maximum *y* to be calculated (it will be ignored if *slice_axis_name* = *y*);
- ***ypoints***: number of points in *y* (it will be ignored if *slice_axis_name* = *y*);
- ***zmin***: minimum *z* to be calculated (it will be ignored if *slice_axis_name* = *z*);
- ***zmax***: maximum *z* to be calculated (it will be ignored if *slice_axis_name* = *z*);
- ***zpoints***: number of points in *z* (it will be ignored if *slice_axis_name* = *z*);
- ***method***: the generation method. A list containing all the methods is available in [[6](#references)]; 
- ***w0***: the gaussian waist for the *method*s *realh_resistant_besselgauss_nonperiodic* and *realh_resistant_besselgauss* (it will be ignored for any other *method*);
- ***file_bbar***: file path in <a href="https://math.nist.gov/MatrixMarket/formats.html">MTX</a> to load the *&#946;&#x035E;*<sub>*p*</sub> parameter responsible to compensate the losses for the method *complexh_resistant_besselgauss* (it will be ignored for any other *method*);
- ***file_sip***: file path to load the SIP *F*(*x*<sub>*p*</sub>,*z*) in <a href="https://math.nist.gov/MatrixMarket/formats.html">MTX</a>;
- ***file_psi***: file path to save the field *Ψ*(*x*,*y*,*z*) in <a href="https://math.nist.gov/MatrixMarket/formats.html">MTX</a>;

An example of usage within <a href="https://www.wolfram.com/mathematica/">Mathematica software</a> is in the zip file *example-with-mathematica.zip* where *scalar-surface-frozen-waves.exe* was compiled at a Windows-x64 architecture.

## References
[1]<a href="https://doi.org/10.1364/JOSAB.412756"> J. O. de Sarro, L. A. Ambrosio, “Surface beams resistant to diffraction and attenuation and structured at the millimeter scale”, 2021.</a>

[2]<a href="https://doi.org/10.1364/JOSAB.36.000638"> L. A. Ambrosio, “Millimeter-structured nondiffracting surface beams”, J. Opt. Soc. Am. B, vol. 36, no. 3, p. 638, Mar. 2019.</a>

[3]<a href="https://doi.org/10.1109/PIERS-Spring46901.2019.9017377"> J. O. de Sarro and L. A. Ambrosio, “Constructing Millimeter-structured Surface Beams from Nondiffracting Zeroth-order Bessel Beams in Lossless Media,” in 2019 PhotonIcs & Electromagnetics Research Symposium - Spring (PIERS-Spring), 2019, pp. 283–288.</a>

[4]<a href="https://doi.org/10.1364/OE.14.001804"> M. Zamboni-Rached, “Diffraction-Attenuation resistant beams in absorbing media,” Opt. Express, vol. 14, no. 5, pp. 1804–1809, Mar. 2006.</a>

[5]<a href="https://doi.org/10.1103/PhysRevA.92.043839"> M. Zamboni-Rached and M. Mojahedi, “Shaping finite-energy diffraction- and attenuation-resistant beams through Bessel-Gauss–beam superposition,” Phys. Rev. A, vol. 92, no. 4, p. 043839, Oct. 2015.</a>

[6] List of the *method*s:
- *complexh_resistant*: a traditional SFW that appears in [[1-4](#references)];
- *complexh_nonresistant*: a non-resistant version of the traditional SFW as described in section 3.A of [[1](#references)];
- *realh_resistant*: method developed in chapter II of [[5](#references)];
- *complexh_resistant_besselgauss*: a modified version of the *complexh_resistant* method where one can control the losses by means of the parameter *&#946;&#x035E;*<sub>*p*</sub>;
- *realh_resistant_besselgauss*: method developed in chapter III of [[5](#references)];
- *realh_resistant_besselgauss_nonperiodic*: method developed in chapter IV of [[5](#references)].;