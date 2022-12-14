#!/bin/bash
cd ~
git clone https://github.com/gjhunt/hspe.git
cd hspe
echo "import(Matrix)" >> lib_hspe/NAMESPACE
sed -i "s/hspe/hspeSparse/g" lib_hspe/DESCRIPTION
sed -i "s/DEoptimR,/DEoptimR,\n    Matrix,/g" lib_hspe/DESCRIPTION
mv lib_hspe hspeSparse
tar -czvf hspeSparse_0.1.tar.gz hspeSparse
