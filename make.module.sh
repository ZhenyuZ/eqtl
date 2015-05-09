#!bin/bash/
rm module.r
touch module.r
cat ./access/*r >> module.r
cat ./annotation/*r >> module.r

