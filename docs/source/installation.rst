Installation
============

Code
------------
#. Clone this repository::

	git clone https://github.com/gait-tech/gaittoolbox.git

#. Run the following command to initialize submodules::

	git submodule update --init --recursive

Libraries
---------
- At least MATLAB 2019
- `MATLAB Aerospace toolbox <https://au.mathworks.com/help/aerotbx/index.html?s_tid=CRUX_lftnav>`_.
- Lie Group SE(3) library (already in repo ```mod-lib/liese3lib```) by `Barfoot and Furgale <http://asrl.utias.utoronto.ca/code/>`_ `doi:10.1109/TRO.2014.2298059 <http://dx.doi.org/10.1109/TRO.2014.2298059>`_
- [Optional] btk toolbox to enable C3D export `Enable C3D Export`_

Enable C3D Export
^^^^^^^^^^^^^^^^^

#. Download the btk toolbox for MATLAB `here <https://code.google.com/archive/p/b-tk/downloads>`_ (e.g., btk-0.3.0_Win7_MatlabR2009b_64bit.zip). 
#. Export to a directory and add the said directory to the MATLAB path. See `tutorial here <https://au.mathworks.com/help/matlab/matlab_env/add-remove-or-reorder-folders-on-the-search-path.html>`_.

.. _Mokka: https://biomechanical-toolkit.github.io/mokka/

After running the code and exporting the C3D, you can visualize the c3d file via the following steps.

#. Download `Mokka`_ and install.
#. Download and load the corresponding **.mvc** file from `gtb Mokka <https://github.com/gait-tech/gaittoolbox/tree/master/mod-lib/Mokka>`_ (see steps in image below). Without loading this config, you will only see the dots but without the colored lines and planes.

	.. image:: fig/mokka-config.png
		:width: 500px
		:align: center
		:alt: Load config to Mokka.
#. See sample output below.

	.. image:: fig/sample-c3d.PNG
		:width: 500px
		:align: center
		:alt: Sample c3d output

Minimal Example
---------------

- See `sample code <https://github.com/gait-tech/gaittoolbox/tree/master/%2Bexamples>`_. 
- Specially recommend checking `runSample02.m <https://github.com/gait-tech/gaittoolbox/blob/master/%2Bexamples/runSample02.m>`_ which contains a simple step by step explanation of how to get body imu measurements with respect the world frame, and then feed it into the CKF filter.

1. For a quick run, go to the root of the gaittoolbox code and run the following command from console. This command will run the `sample code <https://github.com/gait-tech/gaittoolbox/tree/master/%2Bexamples>`_.

.. code-block::

	examples.runSample01


2. To run the experiments in `+paper`, you will run the following command from console

.. code-block::

	papers.<foldername>.runAllNeuRASparse01
	papers.ckf2019.runAllNeuRASparse01