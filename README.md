### SALEc

A implement of SALE algorithm for Crater Formation

This code have been used in Huacheng Li, Zongyu Yue, Yangting Lin, Kaichang Di, Nan Zhang, Jianzhong Liu, Olivine origination in lunar Das crater through three-dimensional numerical simulation, Icarus, 2022 (http://dx.doi.org/10.1016/j.icarus.2022.115333). We recommend authors cite this research in published works.

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
cp ../SALEc.inp .
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
