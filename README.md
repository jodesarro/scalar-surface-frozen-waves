# Surface Frozen Waves Intensity
A C++ routine to evaluate the intensity of scalar zeroth order surface frozen waves (SFWs) [[1-5](#references)].

The ordinary SFW, composed by linear frozen waves (LFWs), posses height and length proportional to the desired superficial intensity pattern (SIP).
The LFWs are equally spaced, with the same parameters L and N but distincts Q<sub>p</sub>.

## Usage
After compiling, the program must be executed with the following parameters passed as arguments, '*nr* *ni* *l0* *file_Qp* *N* *L* *P* *slice_axis_name* *slice_axis_value* *xmin* *xmax* *xpoints* *ymin* *ymax* *ypoints* *zmin* *zmax* *zpoints* *method* *w0* *file_sip* *file_psi*', where:

- ***nr***: real part of the refractive index;
- ***ni***: imaginary part of the refractive index;
- ***l0***: wavelength in vacuum;
- ***file_Qp***: the file path to the *Q<sub>p</sub>* in <a href="https://math.nist.gov/MatrixMarket/formats.html">Matrix Market Array File Format (MTX)</a>;
- ***N***: number related to Bessel beams, resulting in a total of 2*N*+1 Bessel beams;
- ***L***: parameter *L* of the SFWs;
- ***P***: number of LFWs;
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
- ***w0***: the gaussian waist for the *method*s *finite_energy* and *paraxial_apodized* (it will be ignored for any other *method*);
- ***file_sip***: the file path to the SIP *F*(*x*<sub>*p*</sub>,*z*) in <a href="https://math.nist.gov/MatrixMarket/formats.html">MTX</a>;
- ***file_psi***: the file path to save the intensity |*Ψ*(*x*,*y*,*z*)|² in <a href="https://math.nist.gov/MatrixMarket/formats.html">MTX</a>;

An example of usage within <a href="https://www.wolfram.com/mathematica/">Mathematica software</a> is in the zip file *example-with-mathematica.zip* where *surface-frozen-waves.exe* was compiled at a Windows-x64 architecture.

## References
[1] J. O. de Sarro, L. A. Ambrosio, “Surface beams resistant to diffraction and attenuation and structured at the millimeter scale”, 2020.

[2]<a href="https://doi.org/10.1364/JOSAB.36.000638"> L. A. Ambrosio, “Millimeter-structured nondiffracting surface beams”, J. Opt. Soc. Am. B, vol. 36, no. 3, p. 638, Mar. 2019.</a>

[3]<a href="https://doi.org/10.1109/PIERS-Spring46901.2019.9017377"> J. O. de Sarro and L. A. Ambrosio, “Constructing Millimeter-structured Surface Beams from Nondiffracting Zeroth-order Bessel Beams in Lossless Media,” in 2019 PhotonIcs & Electromagnetics Research Symposium - Spring (PIERS-Spring), 2019, pp. 283–288.</a>

[4]<a href="https://doi.org/10.1364/OE.14.001804"> M. Zamboni-Rached, “Diffraction-Attenuation resistant beams in absorbing media,” Opt. Express, vol. 14, no. 5, pp. 1804–1809, Mar. 2006.</a>

[5]<a href="https://doi.org/10.1103/PhysRevA.92.043839"> M. Zamboni-Rached and M. Mojahedi, “Shaping finite-energy diffraction- and attenuation-resistant beams through Bessel-Gauss–beam superposition,” Phys. Rev. A, vol. 92, no. 4, p. 043839, Oct. 2015.</a>

[6] List of the *method*s:
- *traditional*: a traditional SFW that appears in [[1-4](#references)];
- *traditional_nonresistant*: a non-resistant version of the traditional SFW as described in section 3.A of [[1](#references)];
- *modified*: method developed in chapter II of [[5](#references)];
- *paraxial_apodized*: method developed in chapter III of [[5](#references)];
- *finite_energy*: method developed in chapter IV of [[5](#references)].