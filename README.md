### SALEc: 用于撞击坑数值模拟实现的简化任意拉格朗日-欧拉法 [en](README_en.md)
SALEc是基于任意拉格朗日-欧拉法(Hirt et al., 1974; Amsden, 1980)开发的超高速撞击模拟程序。SALEc采用了许多和[iSALE-2D](https://isale-code.github.io/)类似的算法。一些模块的实现受到了iSALE-2D的启发。它可以求解多种物质的数值模型，其中网格内边界建立的算法包括几何法(例如Benson, 2002)和代数法(例如Ubbink and Issa, 1999)，并且使用了MPI并行算法以提高计算效率。SALEc实现了多种本构方程，以及Tillotson状态方程(Tillotson, 1962)和ANEOS状态方程(Thompson and Lauson, 1972)。为了考虑温度和破裂对强度的影响，SALEc采用了与iSALE-2D相同的模型，包括热软化模型(Ohnaka, 1995)和破裂模型(Melosh et al., 1992; Ivanov et al., 1997; Collins et al., 2004)；为了模拟大型撞击坑的垮塌，SALEc还考虑了声波液化模型(Melosh et al., 1979; Melosh and Ivanov, 1999; Collins et al., 2016)，这部分与iSALE-2D也基本相同。

代码已经在Huacheng Li, Zongyu Yue, Yangting Lin, Kaichang Di, Nan Zhang, Jianzhong Liu, Olivine origination in lunar Das crater through three-dimensional numerical simulation, Icarus, 2022 (http://dx.doi.org/10.1016/j.icarus.2022.115333) 中使用。使用SALEc时推荐引用这篇文章。

#### 二维版本[SALEc-2D]([huachengli/SALEc-2D-public](https://github.com/huachengli/SALEc-2D-public))

#### 参考文献
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

1.1 编译

标准的cmake/make编译流程，需要提前安装cmake/make/gcc/mpi工具

```bash
cd salec-download-dir/SALEc
mkdir build
cd build
# change the CMakeLists.txt if needed
cmake ..
make

```

1.2 运行
运行目录的上一层目录中必须要有eos文件夹和参数文件SALEc.inp

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
-np后面的参数必须和$npgx\times npgy\times npgz$相等。

1.3 可视化

可以使用[ParaView](https://www.paraview.org/)或者VisIt打开vtm/vts文件。输出文件太大时也可以使用[SALEcVtsReader](https://github.com/huachengli/SALEcVtsReader)

#### 文档

- 暂时没有，以后可能会加。

#### 作者

- 李华成, huacheng_li@pku.edu.cn
