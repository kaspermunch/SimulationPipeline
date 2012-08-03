
DIR="`pwd``date "+/coalhmm-bpp-%d-%m-%Y"`"

mkdir $DIR
mkdir $DIR/local
mkdir $DIR/local/lib
mkdir $DIR/local/include

export CPATH=$DIR/local/include    # for compilation
export LIBRARY_PATH=$DIR/local/lib    # for compilation
export LD_LIBRARY_PATH=$DIR/local/lib

export C_INCLUDE_PATH=$DIR/local/include
export CPLUS_INCLUDE_PATH=$DIR/local/include

cd $DIR
git clone http://biopp.univ-montp2.fr/git/bpp-core
git clone http://biopp.univ-montp2.fr/git/bpp-seq
git clone http://biopp.univ-montp2.fr/git/bpp-phyl

mkdir ./local
echo `pwd`
cd bpp-core
for package in bpp-core bpp-seq bpp-phyl; do
    cd ../$package
    cmake . -DCMAKE_INSTALL_PREFIX=../local
    make
    make install
#    make apidoc
done

cd .. 
svn co svn://svn.gna.org/svn/coalhmm/trunk coalhmm
cd coalhmm
cmake -DCMAKE_INSTALL_PREFIX=$DIR/local
make
make install


svn co svn://svn.gna.org/svn/bppsuite/trunk bppsuite
cd bppsuite
cmake -DCMAKE_INSTALL_PREFIX=$DIR/local
make
make install

# curl http://download.gna.org/bppsuite/sources/bppsuite-0.6.0.tar.gz -o bppsuite-0.6.0.tar.gz
# tar xfvz bppsuite-0.6.0.tar.gz
# cd bppsuite-0.6.0
# cmake -DCMAKE_INSTALL_PREFIX=../local


cd ..
mkdir package
find ./local -name '*.dylib' -exec cp {} package/. \;
cp ./local/bin/coalhmm package/.

# so when running coalhmm on OSX on grid just include dependencies like this:

#dpu.addDependencies(glob.glob("coalhmm-bpp-17-11-2011/package/*"))


