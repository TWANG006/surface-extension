## Introduction
Different surface extension algorithms for computer controlled optical surfacing (CCOS) are implemented. 

## Implemented algorithms
Surface extension algorithms for arbitrary shape of the original surface error map, including
- [x] Zero extension
- [x] Gaussian extension [3-9]
- [x] Nearest neighbor extension [2, 6, 8]
- [x] Nearest neighbor extension with fall on edge [8]
- [x] C1 smooth extension [5, 6, 8]
- [x] C1 smooth extension withf all on edge [5, 6, 8]
- [x] Gerchberg extension [1, 4, 6]

## Example results
![Surface Extension Results](/images/surface_extension_results.png)

## Usage
**Note:**
The main function is ```[X_ext, Y_ext, Z_ext, ca_range] = Surface_Extension(X, Y, Z, brf_params, Z_tif,method, isFall, fx_range, fy_range)```.
- Inputs:
  - ```X, Y, Z```: coordinate grids and initial surface error map
  - ```brf_params```: Tool Influence Function (TIF) or Beam Removal Function (BRF) parameters, including the Peak Removal Rate (PRR) ```A```, ```sigma_xy```, diameter ```d```, and lateral resolution ```lat_res_brf```
  - ```Z_tif```: Real measured (or simulated) TIF profile
  - ```method```: extension method, including ```zero```, ```gauss```, ```8nn```, ```smooth```, and ```gerchberg```
  - ```isfall```: if fall profile is applied on the extented edge
  - ```fx_range, fyrange```: frequency ranges for Gerchberg's algorithm
- Outputs:
  - ```X_ext, Y_ext, Z_ext```: extended coordinate grids and surface error map
  - ```ca_range```: Range of original surface in the extended surface

## Reference

[1] [Marks, R. J. (1981). Gerchberg’s extrapolation algorithm in two dimensions. Applied optics, 20(10), 1815-1820.](https://doi.org/10.1364/AO.20.001815)

[2] [周林. (2008). 光学镜面离子束修形理论与工艺研究 (Doctoral dissertation, 长沙: 国防科技大学).](http://cdmd.cnki.com.cn/Article/CDMD-90002-2010164869.htm)

[3] [焦长君. (2008). 光学镜面离子束加工材料去除机理与基本工艺研究 (Doctoral dissertation, 长沙: 国防科技大学).](http://cdmd.cnki.com.cn/Article/CDMD-90002-2009213234.htm)

[4] [Liu, Y., Cheng, H., Dong, Z., & Tam, H. Y. (2014). Edge effect of optical surfacing process with different data extension algorithms. Frontiers of Optoelectronics, 7(1), 77-83.](https://link.springer.com/content/pdf/10.1007/s12200-014-0393-7.pdf)

[5] [Zhou, L., Dai, Y., Xie, X., Peng, X., Shi, F., & Li, S. (2016, May). Methods to extend surface error map in dwell time algorithm. In EUSPEN’s 16th International Conference & Exhibition.](https://www.euspen.eu/euspen-knowledge-base/proceedings/)

[6] [Yang, B., Xie, X., Li, F., & Zhou, L. (2017). Edge effect correction using ion beam figuring. Applied Optics, 56(32), 8950-8958.](https://doi.org/10.1364/AO.56.008950)

[7] [Cheng, H. (2016). Independent variables for optical surfacing systems. Springer-Verlag Berlin An.](https://link.springer.com/content/pdf/10.1007/978-3-642-45355-7.pdf)

[8] [唐才学, 颜浩, 罗子健, 张远航, & 温圣林. (2019). 连续位相板磁流变加工中高精度边缘延拓技术. 红外与激光工程, (2019 年 04), 154-160.](https://www.airitilibrary.com/Publication/alDetailedMesh?docid=hwyjggc201904023)

[9] [Zhou, L., Dai, Y., Xie, X., & Li, S. (2016, October). Computations involved in IBF process and introduction of the software IBFCAM. In Advanced Optical Design and Manufacturing Technology and Astronomical Telescopes and Instrumentation (Vol. 10154, p. 101541S). International Society for Optics and Photonics.](https://doi.org/10.1117/12.2247130)
