python3 -m pip install . --no-deps -vv 

mkdir -p "${PREFIX}/bin"

cp -r nanoreceptorModules/* $PREFIX/bin/

