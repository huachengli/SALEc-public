### SALEc
A implement of SALE algorithm for Crater Formation

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

See the document in doc dir

#### Author

- Li Huacheng, huacheng_li@pku.edu.cn

