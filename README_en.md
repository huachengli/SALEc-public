### SALEc: Implement of SALE algorithm for Crater Formation [cn](README.md)

SALEc is a hypervelocity impact simulation program developed based on the Arbitrary Lagrangian-Eulerian (ALE) method (Hirt et al., 1974; Amsden, 1980) and uses many similar algorithms as [iSALE-2D](https://isale-code.github.io/). The code implementation of some modules (e.g., equations of state, strength, acoustic fluidization) in SALEc is largely inspired by iSALE-2D. SALEc solves numerical models involving multiple materials, with algorithms for establishing boundaries within the grid including the geometric method (e.g., Benson, 2002) and the algebraic method (e.g., Ubbink and Issa, 1999). MPI parallel algorithms are utilized to enhance computational efficiency. SALEc incorporates various constitutive equations, as well as the Tillotson equation of state (Tillotson, 1962) and the ANEOS equation of state (Thompson and Lauson, 1972). To account for the effects of temperature and fracture on material strength, SALEc adopts models same as those in iSALE-2D, including the thermal softening model (Ohnaka, 1995) and the fracture model (Melosh et al., 1992; Ivanov et al., 1997; Collins et al., 2004). Additionally, to simulate the collapse of large impact craters, SALEc incorporates the acoustic fluidization model (Melosh et al., 1979; Melosh and Ivanov, 1999; Collins et al., 2016), which is also largely same as iSALE-2D.


This code have been used in Huacheng Li, Zongyu Yue, Yangting Lin, Kaichang Di, Nan Zhang, Jianzhong Liu, Olivine origination in lunar Das crater through three-dimensional numerical simulation, Icarus, 2022 (http://dx.doi.org/10.1016/j.icarus.2022.115333). We recommend authors cite this research in published works.

#### 2d version[huachengli/SALEc-2D-public](https://github.com/huachengli/SALEc-2D-public)

#### References

1. Amsden, A. A., Ruppel, H. M., Hirt, C. W., 1980. SALE: A simplified ALE computer program for fluid flow at all speeds. Los Alamos National Laboratory Technical Report.
2. Benson, D. J., 2002. Volume of fluid interface reconstruction methods for multi-material problems. Applied Mechanics Reviews, 55(2), 151–165.
3. Collins, G. S., Elbeshausen, D., Davison, T. M., Wünnemann, K., Ivanov, B., Melosh, H. J., 2016. iSALE-Dellen manual. Figshare, 136.
4. Collins, G. S., Melosh, H. J., Ivanov, B., 2004. Modeling damage and deformation in impact simulations. Meteoritics & Planetary Sciences, 39(2), 217–231.
5. Hirt, C., Amsden, A., Cook, J., 1974. An arbitrary Lagrangian-Eulerian computing method for all flow speeds. Journal of Computational Physics, 14(3), 227–253.
6. Ivanov, B., Deniem, D., Neukum, G., 1997. Implementation of dynamic strength models into 2D hydrocodes: Applications for atmospheric breakup and impact cratering. International Journal of Impact - Engineering, 20(1–5), 411–430.
7. Melosh, H. J., 1979. Acoustic fluidization - A new geologic process. Journal of Geophysical Research, 84, 7513–7520.
8. Melosh, H. J., Ivanov, B., 1999. Impact crater collapse. Annual Review of Earth and Planetary Sciences, 27(1), 385–415.
9. Melosh, H. J., Ryan, E. V., Asphaug, E., 1992. Dynamic fragmentation in impacts: Hydrocode simulation of laboratory impacts. Journal of Geophysical Research, 97(E9), 14735–14759.
10. Ohnaka, M., 1995. A shear failure strength law of rock in the brittle-plastic transition regime. Geophysical Research Letters, 22(1), 25–28.
11. Thompson, S. L., Lauson, H. S., 1972. Improvements in the Chart D radiation-hydrodynamic code. III: Revised analytic equations of state. Sandia National Laboratory Technical Report, 714.
12. Tillotson, J. H., 1962. Metallic equations of state for hypervelocity impact. General Atomic Technical Report, 3216.
13. Ubbink, O., Issa, R., 1999. A method for capturing sharp fluid interfaces on arbitrary meshes. Journal of Computational Physics, 153(1), 26–50.



#### Get Start

1.1 Compile

```bash
# cmake, make and mpi are required
cd salec-download-dir/SALEc
mkdir build
cd build
# change the CMakeLists.txt if needed
cmake ..
make

```

1.2 Run

```bash
# copy SALEc binary and input setting file to the SALEc working directory(salec-download-dir/SALEc/work/.)
cd ..
mkdir work
cd work
cp ../build/SALEc .
cp ../SALEc_Chicxulub.inp ./SALEc.inp # chicxulub is a Rock-impact benchmark
# link the EoS file to working directory, make data file and run SALEc
ln -s ../eos ./eos
mkdir pdata
cd pdata
mpirun -np 8 ../SALEc
```

User-defined parameters is in SALEc.inp. The number of progresses specified by option "-np" should equal $npgx\times npgy\times npgz$ in SALEc.inp.

1.3 Visualize

Using [ParaView](https://www.paraview.org/) to open the *.vtm file in pdata.

#### More Information

- document will be added in future

#### Author

- Li Huacheng, huacheng_li@pku.edu.cn
