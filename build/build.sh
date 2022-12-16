#!/bin/sh
set -e

$PYTHON setup.py install

# mkdir -p "${PREFIX}/bin"
# mkdir -p "${PREFIX}/bin/modules"
# mkdir -p "${PREFIX}/databases"

# cp -r bin/* "${PREFIX}/bin/"
# cp -r bin/modules/* "${PREFIX}/bin/modules"
# cp -r databases/* "${PREFIX}/databases/"
