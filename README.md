# QuadrantPackage
Quadrant Detection Package, as described in M. A. Taylor and W. P. Bowen, “A computational tool to characterize particle tracking measurements in optical tweezers,” J. Opt. 15, 085701 (2013). 

version 1.1

Contents of this README file:

1. Warning
2. License
3. Package contents
4. Getting started
5. References

Contact: m.taylor@sbs.uq.edu.au

============
1. Warning
============

This code is a work in progress and may contain bugs. I do not guarantee the 
accuracy of all results. The expansions used can be susceptible to numerical 
errors due to the addition of very many numbers, particularly for very large 
particle size. Please let me know of any bugs you find. 

Any version updates will be posted at
https://github.com/michael-a-taylor/QuadrantPackage

Contact me if you require help getting this to work.



THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE. 

==========
2. License
==========

Copyright 2013-2016 The University of Queensland.

This package and its component files are copyright 2016 by The University of 
Queensland. Not-for-profit re-distribution of the package is permitted. The 
package may be used free-of-charge for research, teaching, or personal use. 

This package includes code from the Optical Tweezers Toolbox version 1.3, 
Copyright 2007-2014 The University of Queensland, reused with permission. 
The Optical Tweezers Toolbox is hosted at: 
http://www.physics.uq.edu.au/people/nieminen/software.html
 
If results obtained using this package are published, both the package and the 
Optical Tweezers Toolbox should be appropriately referenced.


===================
3. Package contents
===================
The code was written as an add-on to the Optical Tweezers Toolbox. Most 
of the files included in the package are copied unmodified from the Optical 
Tweezers Toolbox version 1.3. The Quadrant Detection Package itself only 
comprises the following files: 

- Quadrant_measurement.m: This calculates quadrant detection signal and shot-noise limit to displacement precision for homogeneous spheres. This is the main file to look at in this package.
- Quadrant_measurement_layered.m: This calculates the quadrant detection signal for a layered sphere.
- Quadrant_measurement_cube.m: This calculates the quadrant detection signal for a cube.
-farfield_matrix.m: This is the original file as presented in the paper. It calculates matrices that convert between the E(theta,phi) basis and vector spherical harmonics.
-farfield_matrix2.m: This is a modified version of farfield_matrix.m that is computationally far more efficient.

Additionally, I have modified the file spharm.m to improve numerical stability.

The biggest change since the original release with the paper is the increased speed provided with farfield_matrix2.m. Additionally, there have been some minor bug fixes


==================
4. Getting started
==================


(a) Firstly, read the paper. This will outline the basic principle of the algorithm, and how it is implemented in this software package.

(b) Install the package. Unzip all of the files into a directory. Either work in that directory, or add it to your MATLAB path.

(c) Play around with the code. The only file that needs attention is Quadrant_measurement.m. 


(d) All calculations work in length units of the medium wavelength, but I decided to specify user inputs in meters.
    

=============
5. References
=============

Quadrant Detection Package:
M. A. Taylor and W. P. Bowen, "A computational tool to characterize particle 
tracking measurements in optical tweezers," J. Opt. 15, 085701 (2013). 
https://github.com/michael-a-taylor/QuadrantPackage


The Optical Tweezers Toolbox package:
T. A. Nieminen, V. L. Y. Loke, A. B. Stilgoe, G. Knoener,
A. M. Branczyk, N. R. Heckenberg, H. Rubinsztein-Dunlop,
"Optical tweezers computational toolbox",
Journal of Optics A 9, S196-S203 (2007)

T. A. Nieminen, A. B. Stilgoe, V. L. Y. Loke, N. du Preez-Wilkinson,
 A. A. M. Bui, Y. Cao, Y. Hu, G. Knoener, A. M. Branczyk,
"Optical tweezers toolbox 1.3",
http://www.physics.uq.edu.au/people/nieminen/software.html
