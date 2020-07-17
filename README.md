# Surface Frozen Waves Intensity
A C++ routine to evaluate the intensity of Surface Frozen Waves (SFWs) [[1](https://github.com/jodesarro/sfwi#references)].
The code is a structured routine that may need improvement.
It is not a "Clean Code", but it works!

## Usage
After compiling, the program must be executed with the following parameters passed as arguments, 'nr ni l0 Q N L P cut_axis cut_value xmin xmax xpoints ymin ymax ypoints zmin zmax zpoints is_nonresistant file_in file_out', where:

- **nr**: real part of the refractive index;
- **ni**: imaginary part of the refractive index;
- **l0**: wavelength in vacuum;
- **Q**: parameter Q of the SFWs;
- **N**: number related to Bessel beams, resulting in a total of 2N+1 Bessel beams;
- **L**: parameter L of the SFWs;
- **P**: number of linear FWs;
- **cut_axis**: axis x, y or z where the sliceplot will be perpendicular to it;
- **cut_value**: value to slice the cut_axis;
- **xmin**: minimum x to be calculated (it will be ignored if cut_axis = x);
- **xmax**: maximum x to be calculated (it will be ignored if cut_axis = x);
- **xpoints**: number of points in x (it will be ignored if cut_axis = x);
- **ymin**: minimum y to be calculated (it will be ignored if cut_axis = y);
- **ymax**: maximum y to be calculated (it will be ignored if cut_axis = y);
- **ypoints**: number of points in y (it will be ignored if cut_axis = y);
- **zmin**: minimum z to be calculated (it will be ignored if cut_axis = z);
- **zmax**: maximum z to be calculated (it will be ignored if cut_axis = z);
- **zpoints**: number of points in z (it will be ignored if cut_axis = z);
- **is_nonresistant**: 0 for resistant SFWs and 1 for nonresistant SFWs;
- **file_in**: the file path to the F function in <a href="https://math.nist.gov/MatrixMarket/formats.html">Matrix Market Array File Format (MTX)</a>;
- **file_out**: the file path to save the intensities in <a href="https://math.nist.gov/MatrixMarket/formats.html">MTX</a>;

An example of usage within <a href="https://www.wolfram.com/mathematica/">Mathematica software</a> is in the folder [example-with-mathematica](example-with-mathematica) where sfwi.exe was compiled at a Windows-x64 architecture.

## References
[[1](https://github.com/jodesarro/sfwi#surface-frozen-waves-intensity)] <a href="https://doi.org/10.1364/JOSAB.36.000638">L. A. Ambrosio, “Millimeter-structured nondiffracting surface beams,” J. Opt. Soc. Am. B, vol. 36, no. 3, p. 638, Mar. 2019.</a>
