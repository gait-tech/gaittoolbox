Installation
============

Core
------------

1. Clone this repository 

.. code-block::

	git clone https://github.com/gait-tech/gaittoolbox.git

2. Run the following command to initialize submodules

.. code-block::

	git submodule update --init --recursive

	
Enable C3D Export
-----------------
#. Download the btk toolbox for MATLAB `here <https://code.google.com/archive/p/b-tk/downloads>`_ (e.g., btk-0.3.0_Win7_MatlabR2009b_64bit.zip). 
#. Export to a directory and add the said directory to the MATLAB path. See `tutorial here <https://au.mathworks.com/help/matlab/matlab_env/add-remove-or-reorder-folders-on-the-search-path.html>`_.

Minimal Example
---------------

1. For a quick run, go to the root of the gaittoolbox code and run the following command from console. This command will run the `sample code <https://github.com/gait-tech/gaittoolbox/tree/master/%2Bexamples>`_.

.. code-block::

	examples.runSample01


2. To run the experiments in `+paper`, go to the root of the gaittoolbox code and run the following command from console.

.. code-block::

	papers.<foldername>.runAllNeuRASparse01
	papers.ckf2019.runAllNeuRASparse01
