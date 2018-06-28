## How to start
1. Install MATLAB 2018a (there are basic functions that uses the newest library)
1. Run runAllTCD.m, runAllXsens.m, OR runAllNeuRA.m (not yet ready)
1. To run UKF based algorithms, refer to UKF_validate_cukf_v9.m or similar files

## Total Capture Dataset
1. Download from [link](https://unsw-my.sharepoint.com/:u:/g/personal/z5151460_ad_unsw_edu_au/EeKUGznFC3tAjUneixJmh64B93ozvw5uKnYt9gl2KecDxw?e=EjeVeP). Please do not forget to cite the reference paper [details](http://cvssp.org/data/totalcapture/)
1. Extract zip file to the code root folder. It should produce a folder named 'totalcapture' with subdirectories 'gyroMag', 'imu', and 'vicon'.

## Reference:
- BVHBody's loadbvh under tcdlib is taken from https://github.com/wspr/bvh-matlab