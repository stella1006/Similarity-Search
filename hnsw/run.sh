vecsize=17770
vecdim=300
dataset="netflix"

#vecsize=2340373
#vecdim=150
#dataset="imagenet"

prefix="/data/jinfeng/project/github/gqr/data/${dataset}"
lshboxPath="${prefix}/${dataset}_product_groundtruth.lshbox"
basePath="${prefix}/${dataset}_base.fvecs"
queryPath="${prefix}/${dataset}_query.fvecs"

cmake ./ 
make  2>&1 | tee log.txt
time ./main $vecsize $vecdim $lshboxPath $basePath $queryPath
