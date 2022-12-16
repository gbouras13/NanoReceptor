mkdir -p "${PREFIX}/bin/nanoreceptorModules"

cp -r nanoreceptorModules/* $PREFIX/bin/nanoreceptorModules/

python -m pip install . --no-deps -vv 