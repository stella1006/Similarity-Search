
topk=10
construction_k=20

# 0 -> construction ; 1 -> search
mode=$1
dataset=$2
if [ "$2" = "netflix" ]; then
    vecsize=17770
    vecdim=300
elif [ "$2" = "yahoomusic" ]; then
    vecsize=136736
    vecdim=300
elif [ "$2" = "imagenet" ]; then
    vecsize=2340373
    vecdim=150
fi


prefix="/data/liujie/gqr/data/${dataset}"
lshboxPath="${prefix}/${dataset}_top${topk}_product_groundtruth.lshbox"
basePath="${prefix}/${dataset}_base.fvecs"
queryPath="${prefix}/${dataset}_query.fvecs"

#cmake ./ 
make  2>&1 | tee log.txt
if [ $1 -eq 0 ]; then
    # construction
    time ./exknn $mode $dataset $vecsize $vecdim $construction_k $basePath 
else
    # search
    time ./exknn $mode $dataset $vecsize $vecdim $construction_k $basePath $queryPath $lshboxPath $topk
fi
