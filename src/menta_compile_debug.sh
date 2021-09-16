mkdir scsm
cp ../o/*/*.o scsm/
rm scsm/wwm_*.o

export LD_LIBRARY_PATH=/APPLICATIONS/hdf5/g_hdf5_1.8.16/lib
export MPI_PATH=/ADAPTATION/mentalo/usr/mpich_gnu/mpich-3.2.1/
export NETCDF_FORTRAN_LINK="/APPLICATIONS/netcdf/g_netcdf_f_4.4.3/lib/libnetcdff.a /APPLICATIONS/netcdf/g_netcdf_c_4.4.0/lib/libnetcdf.so"
export NETCDF_INCDIR=/APPLICATIONS/netcdf/g_netcdf_c_4.4.0/include/
export NETCDFF_INCDIR=/APPLICATIONS/netcdf/g_netcdf_f_4.4.3/include/
export HDF_INCDIR=/APPLICATIONS/hdf5/g_hdf5_1.8.16/include/

export METIS_PATH=../ParMetis-4.0.3/
export PDLIB_PATH=../../../pdlib/
ln -sf makefile.menta_debug makefile
make clean
make

