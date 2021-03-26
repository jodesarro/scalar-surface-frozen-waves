# Surface frozen waves
A C++ routine to evaluate the fields of zeroth order surface frozen waves (SFWs) [[1-5](#references)].

The ordinary SFW is composed by linear frozen waves (LFWs) equally spaced, with the same parameter *N* but distincts *Q*<sub>*p*</sub> and *L*<sub>*p*</sub> and is intended to recreate a pre-chosen surface intensity pattern (SIP).

## Usage
Before compiling, one may change the following variables inside ```int main()``` function as desired:

- ```string pm.method```: generation method. A list containing all the methods is available in [[6](#references)];
- ```string pm.polarization```: polarization of the electromagnetic field. A list containing the polarizations is available in [[7](#references)];
- ```double pm.nr```: real part of the complex refractive index;
- ```double pm.ni```: imaginary part of the complex refractive index;
- ```double pm.l0```: wavelength in vacuum;
- ```double pm.N```: number related to Bessel beams, resulting in a total of 2*N*+1 Bessel beams;
- ```double pm.P```: number of LFWs;
- ```double * L```: parameter *L*<sub>*p*</sub>;
- ```double * Q```: parameter *Q*<sub>*p*</sub>;
- ```double * w```: parameter *w*<sub>*p*</sub> related to Gaussian waist for each LFW for ```pm.method=="finite_energy"``` (it will be ignored for any other method);
- ```double pm.dx0```: separation in *x* between two consecutives LFWs;
- ```double * pm.x0```: displacement in *x* of the LFWs;
- ```double pm.xmin```: minimum *x* to be calculated;
- ```double pm.xmax```: maximum *x* to be calculated;
- ```int pm.xpoints```: number of points in *x*;
- ```double pm.ymin```: minimum *y* to be calculated;
- ```double pm.ymax```: maximum *y* to be calculated;
- ```int pm.ypoints```: number of points in *y*;
- ```double pm.zmin```: minimum *z* to be calculated;
- ```double pm.zmax```: maximum *z* to be calculated;
- ```int pm.zpoints```: number of points in *z*;
- ```bool pm.is_async```: determines whether to use parallel calculation. If ```true```, the calculation speeds up with the price of slowing down the computer;

Before running the compiled, a file named ```sip.mtx``` containing the SIP *F*(*x*<sub>*p*</sub>,*z*) in <a href="https://math.nist.gov/MatrixMarket/formats.html">Matrix Market Array File Format (MTX)</a> must be in the same folder as the executable.

Once the calculation is done, the program will export the components of the electromagnetic vector field ***E***(*x*,*y*,*z*) as ```Ex.m```, ```Ey.m``` and```Ez.m``` in <a href="https://reference.wolfram.com/language/ref/format/WL.html">WL file type</a> (data in plain ASCII text format) compatible also with Mathematica.

A commented example of usage with <a href="https://www.wolfram.com/mathematica/">Mathematica software</a> is in the file ```example-with-mathematica.zip```.

## References
[1]<a href="https://doi.org/10.1364/JOSAB.412756"> J. O. de Sarro, L. A. Ambrosio, “Surface beams resistant to diffraction and attenuation and structured at the millimeter scale”, J. Opt. Soc. Am. B 38, 677-684, 2021.</a>

[2]<a href="https://doi.org/10.1364/JOSAB.36.000638"> L. A. Ambrosio, “Millimeter-structured nondiffracting surface beams”, J. Opt. Soc. Am. B, vol. 36, no. 3, p. 638, Mar. 2019.</a>

[3]<a href="https://doi.org/10.1109/PIERS-Spring46901.2019.9017377"> J. O. de Sarro and L. A. Ambrosio, “Constructing Millimeter-structured Surface Beams from Nondiffracting Zeroth-order Bessel Beams in Lossless Media,” in 2019 PhotonIcs & Electromagnetics Research Symposium - Spring (PIERS-Spring), 2019, pp. 283–288.</a>

[4]<a href="https://doi.org/10.1364/OE.14.001804"> M. Zamboni-Rached, “Diffraction-Attenuation resistant beams in absorbing media,” Opt. Express, vol. 14, no. 5, pp. 1804–1809, Mar. 2006.</a>

[5]<a href="https://doi.org/10.1103/PhysRevA.92.043839"> M. Zamboni-Rached and M. Mojahedi, “Shaping finite-energy diffraction- and attenuation-resistant beams through Bessel-Gauss–beam superposition,” Phys. Rev. A, vol. 92, no. 4, p. 043839, Oct. 2015.</a>

[6] List of the ```pm.method```'s available:
- ```"resistant"```: traditional SFW that appears in [[1-4](#references)];
- ```"nonresistant"```: non-resistant version of the traditional SFW as described in [[3](#references)] and in section 3.A of [[1](#references)];
- ```"realh_resistant"```: method developed in chapter II of [[5](#references)];
- ```"finite_energy"```: method developed in chapter IV of [[5](#references)];

[7] List of the ```pm.polarization```'s available:
- ```"scalar"```: scalar field *ψ*(*x*,*y*,*z*) that will be exported in ```Ex.m```;
- ```"linear"```: linear polarized field in *x* direction;
- ```"linear_crossed"```: linear cross-polarized field alternating between *x* and *y* directions;