# Gait Reconstruction Toolbox

Refer to the [documentation](https://gait-tech.github.io/gaittoolbox/) for more details.

## Installation

1. Clone this repository 
	```
	git clone https://github.com/gait-tech/gaittoolbox.git
	```
2. Run the following command to initialize submodules
	```
	git submodule update --init --recursive
	```
3. *OPTIONAL but required library to generate c3d files.* Download the btk toolbox for MATLAB [here](https://code.google.com/archive/p/b-tk/downloads) (e.g., btk-0.3.0_Win7_MatlabR2009b_64bit.zip). Export to a directory and add the said directory to the MATLAB path [tutorial here](https://au.mathworks.com/help/matlab/matlab_env/add-remove-or-reorder-folders-on-the-search-path.html).

## Update
1. Pull this repository and submodules
	```
	git pull
	git submodule update --remote --merge
	```

## Documentation and Tutorial

* [Documentation](https://gait-tech.github.io/gaittoolbox/)
* See `+examples` for sample codes

## How to cite the Gait Toolbox

When citing the Gait Toolbox, please cite the original paper where an algorithm was first reported, so that contributors of new algorithms get their fair share of citations. For example, the algorithm *Constrainted Kalman filter for reduced sensor configaration* can be cited as follows:

	[1] Luke Sy, Michael Raitor, Michael Del Rosario, Heba Khamis, Lauren Kark, Nigel H. Lovell, Stephen J. Redmond (2019). Estimating Lower Limb Kinematics using a Reduced Wearable Sensor Count. Pre-print.

## Disclaimer
*The software provided is distributed under the GNU GPLv3 or later. However, this software is designed for scientific research and as such may contain algorithms that are associated with patents in Australia and abroad. If the user so chooses to use the software provided for commercial endeavors then it is solely the user’s responsibility to license any patents that may exist and respond in full to any legal actions taken by the patent holder.*

