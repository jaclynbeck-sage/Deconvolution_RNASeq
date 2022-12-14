#!/bin/bash
cd ~
git clone https://github.com/gjhunt/dtangle.git
cd dtangle
echo "import(Matrix)" >> lib_dtangle/NAMESPACE
sed -i "s/dtangle/dtangleSparse/g" lib_dtangle/DESCRIPTION
sed -i "s/DEoptimR,/DEoptimR,\n    Matrix,/g" lib_dtangle/DESCRIPTION
mv lib_dtangle dtangleSparse
tar -czvf dtangleSparse_2.0.9.tar.gz dtangleSparse
