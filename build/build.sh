python3 -m pip install . --no-deps -vv 

mkdir -p "${PREFIX}/bin"
mkdir -p "${PREFIX}/bin/nanoreceptorModules"

cp -r nanoreceptor.py $PREFIX/bin/
cp -r nanoreceptorModules/* $PREFIX/bin/nanoreceptorModules/

