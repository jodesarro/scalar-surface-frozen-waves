# Surface Frozen Waves Intensity
A C++ routine to evaluate the intensity of scalar zeroth order Surface Frozen Waves (SFWs) [[1-5](#references)].

## Usage
After compiling, the program must be executed with the following parameters passed as arguments, '*nr* *ni* *l0* *Q* *N* *L* *P* *slice_axis_name* *slice_axis_value* *xmin* *xmax* *xpoints* *ymin* *ymax* *ypoints* *zmin* *zmax* *zpoints* *method* *w0* *file_in* *file_out*', where:

- ***nr***: real part of the refractive index;
- ***ni***: imaginary part of the refractive index;
- ***l0***: wavelength in vacuum;
- ***Q***: parameter *Q* of the SFWs;
- ***N***: number related to Bessel beams, resulting in a total of 2*N*+1 Bessel beams;
- ***L***: parameter *L* of the SFWs;
- ***P***: number of linear FWs;
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
- ***file_in***: the file path to the *F* function in <a href="https://math.nist.gov/MatrixMarket/formats.html">Matrix Market Array File Format (MTX)</a>;
- ***file_out***: the file path to save the intensities in <a href="https://math.nist.gov/MatrixMarket/formats.html">MTX</a>;

An example of usage within <a href="https://www.wolfram.com/mathematica/">Mathematica software</a> is in the folder [example-with-mathematica](example-with-mathematica) where *surface-frozen-waves.exe* was compiled at a Windows-x64 architecture.

## References
[1]

[2]<a href="https://doi.org/10.1364/JOSAB.36.000638"> L. A. Ambrosio, “Millimeter-structured nondiffracting surface beams”, J. Opt. Soc. Am. B, vol. 36, no. 3, p. 638, Mar. 2019.</a>

[3]<a href="https://doi.org/10.1109/PIERS-Spring46901.2019.9017377"> J. O. de Sarro and L. A. Ambrosio, “Constructing Millimeter-structured Surface Beams from Nondiffracting Zeroth-order Bessel Beams in Lossless Media,” in 2019 PhotonIcs & Electromagnetics Research Symposium - Spring (PIERS-Spring), 2019, pp. 283–288.</a>

[4]<a href="https://doi.org/10.1364/OE.14.001804"> M. Zamboni-Rached, “Diffraction-Attenuation resistant beams in absorbing media,” Opt. Express, vol. 14, no. 5, pp. 1804–1809, Mar. 2006.</a>

[5]<a href="https://doi.org/10.1103/PhysRevA.92.043839"> M. Zamboni-Rached and M. Mojahedi, “Shaping finite-energy diffraction- and attenuation-resistant beams through Bessel-Gauss–beam superposition,” Phys. Rev. A, vol. 92, no. 4, p. 043839, Oct. 2015.</a>

[6] List of the *method*s:
- *traditional*: a traditional SFW that appears in [[1-4](#references)];
- *traditional_nonresistant*: a non-resistant version of the SFW described in section 3.A of [[1](#references)];
- *modified*: method developed in chapter II of [[5](#references)];
- *paraxial_apodized*: method developed in chapter III of [[5](#references)];
- *finite_energy*: method developed in chapter IV of [[5](#references)].
