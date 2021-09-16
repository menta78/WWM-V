mkdir scsm
cp ../o/*/*.o scsm/
rm scsm/wwm_*.o

export NETCDF_FORTRAN_LINK="/usr/lib/x86_64-linux-gnu/libnetcdff.a /usr/lib/x86_64-linux-gnu/libnetcdf.so"
export NETCDF_INCDIR=/usr/include/
export METIS_PATH=../ParMetis-4.0.3/
export PDLIB_PATH=../../../pdlib/
ln -sf makefile.menta_pcmenta makefile
make clean
make

