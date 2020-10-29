Data Set
============

Neura Sparse Dataset
--------------------

Setup
^^^^^
The experiment had nine healthy subjects (7 men and 2 women, weight :math:`63.0 \pm 6.8` kg, height :math:`1.70 \pm 0.06` m, age :math:`24.6 \pm 3.9` years old), with no known gait or lower body biomechanical abnormalities.
Each subject was compared to two benchmark systems, namely the Vicon and Xsens systems. 

- The Vicon Vantage system consisted of eight cameras covering approximately :math:`4 \times 4 m^2` capture area with millimetre accuracy. Vicon data were captured at 100 Hz and processed using Nexus 2.7 software.
- The Xsens Awinda system consisted of seven MTx units (IMUs). Xsens data were captured at 100 Hz using MT Manager 4.8 and processed using MVN Studio 4.4 software.

The Vicon and Xsens recordings were synchronized by having the Xsens Awinda station send a trigger pulse to the Vicon system at the start and stop event of each recording.
Each subject had reflective Vicon markers placed according to the Helen-Hayes 16 marker set, 
seven MTx units attached to the pelvis, thighs, shanks, and feet according to standard Xsens sensor placement, and two MTx units attached near the ankles.
The MTx units were all factory calibrated. 

Refer to the following papers for examples on how the dataset can be used:
:cite:`sy2019ckf`
:cite:`sy2019lgcekf`.
    
Movements
^^^^^^^^^
Each subject performed the movements listed in the table below twice (i.e., two trials). 
The subjects stood still before and after each trial for ten seconds.
The experiment was approved by the Human Research Ethics Board of the University of New South Wales (UNSW) with approval number HC180413.

+-----------------+---------------------------------------+--------------+
| Movement        | Description                           | Duration (s) |
+-----------------+---------------------------------------+--------------+
| Static          | Stand still                           |  ~10         |
+-----------------+---------------------------------------+--------------+
| Walk            | Walk straight and back                |  ~30         |
+-----------------+---------------------------------------+--------------+
| Figure of eight | Walk in figures of eight              |  ~60         |
+-----------------+---------------------------------------+--------------+
| Zig-zag         | Walk zigzag                           |  ~60         |
+-----------------+---------------------------------------+--------------+
| 5-minute walk   | Undirected walk, side step, and stand | ~300         |
+-----------------+---------------------------------------+--------------+
| Speedskater     | Speedskater on the spot               |  ~30         |
+-----------------+---------------------------------------+--------------+
| Jog             | Jog straight and return               |  ~30         |
+-----------------+---------------------------------------+--------------+
| Jumping jacks   | Jumping jacks on the spot             |  ~30         |
+-----------------+---------------------------------------+--------------+
| High knee       | High knee jog straight and return     |  ~30    	 |
+-----------------+---------------------------------------+--------------+

Installation
^^^^^^^^^^^^
#. Download `neura-sparse01.zip` from `Harvard Dataverse <https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/9QDD5J>`_. Refer to :ref:`desc-base` for details.
#. Extract to `data/neura-sparse01`
#. If you want the raw recordings (not time aligned, and no calibration done), download `neura-sparse01-raw.zip` from the link above. Refer to :ref:`desc-raw` for details.

.. _desc-base:

Description (base)
^^^^^^^^^^^^^^^^^^
<Subject ID>-Trial-<Movement Type>.<Extension>

Note that imu, vicon, and step-detect are all aligned (time-wise), and are with respect the world frame.

imu
"""

File name is <Subject ID>-Trial-<Movement Type>-<BodyPart>.csv.


+--------------+-------------------------------------+
| Body Segment | Sensor Location                     |
+--------------+-------------------------------------+
| Pelvis       | Pelvis  backside                    |
+--------------+-------------------------------------+
| L_UpLeg      | Left thigh                          |
+--------------+-------------------------------------+
| R_UpLeg      | Right thigh                         |
+--------------+-------------------------------------+
| L_LowLeg     | Left ankle (slightly above)         |
+--------------+-------------------------------------+
| R_LowLeg     | Right ankle (slightly above)        |
+--------------+-------------------------------------+
| L_LowLeg2    | Left shanks                         |
+--------------+-------------------------------------+
| R_LowLeg2    | Right shanks                        |
+--------------+-------------------------------------+
| L_Foot       | Left foot  (w/ shoes)               |
+--------------+-------------------------------------+
| R_Foot       | Right foot (w/ shoes)               |
+--------------+-------------------------------------+

Column Description:

+----------------+-------+------------------------------------------------+
| Column         | Type  | Description                                    |
+----------------+-------+------------------------------------------------+
| Acc_X/Y/Z      | Float | Measured acceleration in sensor frame          |
+----------------+-------+------------------------------------------------+
| Gyr_X/Y/Z      | Float | Measured angular velocity in sensor frame      |
+----------------+-------+------------------------------------------------+
| Mag_X/Y/Z      | Float | Measured magnetic field in sensor frame        |
+----------------+-------+------------------------------------------------+
| ori_q0-4       | Float | Orientation of sensor in world frame (w,x,y,z) |
+----------------+-------+------------------------------------------------+

Example::

	PacketCounter,SampleTimeFine,Acc_X,Acc_Y,Acc_Z,Gyr_X,Gyr_Y,Gyr_Z,Mag_X,Mag_Y,Mag_Z,Quat_q0,Quat_q1,Quat_q2,Quat_q3
	27626,,8.696773,0.341363,4.715198,0.015417,0.009366,-0.011596,0.886719,0.420898,0.468506,0.651925,-0.318177,-0.399341,-0.560611
	...
	
vicon
"""""
CSV export of ViconBody. Contains marker and joint centre positions of the movement trial.

xsens
"""""
Pose reconstruction by the Xsens system. Specifically, bvh files of the movement trial generated by the MVN Studio 4.4 software. Note that its index are not aligned with vicon and imu (haven't edited them yet). Data loader code for imu should return the corresponding index that must be used. This data is for comparing with the xsens MVN black box output.

step-detect
"""""""""""
Indicate if step is detected. Reviewed manually.

+--------+---------+-------------------------------------------+
| Column | Type    | Description                               |
+--------+---------+-------------------------------------------+
| stepL  | Boolean | 1 if left foot step is detected, 0 if not |
+--------+---------+-------------------------------------------+
| stepR  | Boolean | 1 if left foot step is detected, 0 if not |
+--------+---------+-------------------------------------------+

Example::

	stepL,stepR
	1,1
	...

calib
"""""

- <Subject ID>-Trial-<Movement>-Calib-SensorYawFixWorldFrame.txt: Contains yaw offset calibration for pelvis, ankle, and foot IMUs. Only important file in the folder.
- <Subject ID>-Calib-V2W-Compass.mat and <Subject ID>-Calib-V2W-Pendumum.mat: used to calculate Vicon to World rotation matrix. Only used in raw processing.
- <Subject ID>-Calib-W2V.txt: Description to follow. Did not use in the dataset.

.. _desc-raw:

Description (raw)
^^^^^^^^^^^^^^^^^

rawvicon
""""""""
CSV export from vicon but converted to .mat file. Contains marker and joint centre positions of the movement trial.

rawimu
""""""
File name is <Subject ID>-Trial-<Movement Type>-<Sensor ID>.txt.

Sensor ID to body segment table:

+--------------+-----------+
| Body Segment | Sensor ID |
+--------------+-----------+
| Pelvis       | 00B40B91  |
+--------------+-----------+
| L_UpLeg      | 00B40C45  |
+--------------+-----------+
| R_UpLeg      | 00B40C3C  |
+--------------+-----------+
| L_LowLeg     | 00B40C44  |
+--------------+-----------+
| R_LowLeg     | 00B40C47  |
+--------------+-----------+
| L_LowLeg2    | 00B40BA5  |
+--------------+-----------+
| R_LowLeg2    | 00B40C35  |
+--------------+-----------+
| L_Foot       | 00B40C55  |
+--------------+-----------+
| R_Foot       | 00B40C48  |
+--------------+-----------+

Column Description:

+----------------+-------+-------------------------------------------+
| Column         | Type  | Description                               |
+----------------+-------+-------------------------------------------+
| PacketCounter  | Int   | Packet number                             |
+----------------+-------+-------------------------------------------+
| SampleTimeFine | Float | Time of recording                         |
+----------------+-------+-------------------------------------------+
| Acc_X/Y/Z      | Float | Measured acceleration in sensor frame     |
+----------------+-------+-------------------------------------------+
| Gyr_X/Y/Z      | Float | Measured angular velocity in sensor frame |
+----------------+-------+-------------------------------------------+
| Mag_X/Y/Z      | Float | Measured magnetic field in sensor frame   |
+----------------+-------+-------------------------------------------+
| Quat_q0-4      | Float | Orientation of sensor in world frame      |
+----------------+-------+-------------------------------------------+

Example::

	PacketCounter,SampleTimeFine,Acc_X,Acc_Y,Acc_Z,Gyr_X,Gyr_Y,Gyr_Z,Mag_X,Mag_Y,Mag_Z,Quat_q0,Quat_q1,Quat_q2,Quat_q3
	27626,,8.696773,0.341363,4.715198,0.015417,0.009366,-0.011596,0.886719,0.420898,0.468506,0.651925,-0.318177,-0.399341,-0.560611
	...
	
rawstep-detect
""""""""""""""
Indicate if step is detected for the whole raw trial. Reviewed manually.