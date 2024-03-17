sudo apt install autoconf

wget --progress=dot:mega https://ftp.gnu.org/gnu/gmp/gmp-6.3.0.tar.xz ; 
tar xf gmp-6.3.0.tar.xz ; cd gmp-6.3.0 ; 
./configure --enable-cxx --enable-fat --prefix=$(pwd)/../usr/local  --build=westmere-pc-linux-gnu ; m
ake -j ; 
make install ; 
cd .. ;
rm gmp-6.3.0.tar.xz ;

wget --progress=dot:mega https://github.com/libexpat/libexpat/archive/R_2_2_4.tar.gz ; 
tar xzf R_2_2_4.tar.gz ; 
cd libexpat-R_2_2_4/expat/ ; 
./buildconf.sh ; 
./configure --prefix=$(pwd)/../../usr/local --without-xmlwf ; 
make -j ; 
make install ; 
cd ../.. ;
rm R_2_2_4.tar.gz ;

wget --progress=dot:mega https://github.com/lip6/libDDD/raw/gh-pages/linux.tgz ; 
tar xzf linux.tgz ; 
rm linux.tgz ;

wget --progress=dot:mega https://yanntm.github.io/Spot-BinaryBuilds/spot_linux.tar.gz ; 
tar xzf spot_linux.tar.gz ;
rm spot_linux.tar.gz ;

cd Petri ; 
autoreconf -vfi ;
./configure --prefix=$PWD/../usr/local/  --with-libexpat=$PWD/../usr/local/ --with-libspot=$(pwd)/../usr/local/ --with-libddd=$PWD/../usr/local/ --with-gmp=$PWD/../usr/local/ --with-antlrc=$PWD/../usr/local/   CPPFLAGS="-I$(pwd)/../usr/local/include -DNDEBUG" LDFLAGS="-L$(pwd)/../usr/local/lib" || cat config.log ;
make -j 4 ;
