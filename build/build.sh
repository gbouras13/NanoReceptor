mkdir -p "${PREFIX}/bin/nanoreceptorModules"

cp -r nanoreceptorModules/* $PREFIX/bin/nanoreceptorModules/

python3 -m pip install . --no-deps -vv 