#include <iostream>
#include <vector>
#include <utility>
#include <cstring>
#include <queue>
#include <ctime>
#include <memory>
#include <algorithm>
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */
#include <omp.h>
#include <chrono>
#include <data.h>
#include <faiss/IndexFlat.h>
#include <faiss/IndexIVFFlat.h>
#include <faiss/IndexPQ.h>
#include <faiss/index_io.h>
using namespace std;

typedef Dataset<float> FloatDataset;

struct IncreaseCmp
{
    bool operator()(const pair<int,float> &a, const pair<int,float> &b)
    {
        return a.second > b.second;
    }
};

struct DecreaseCmp
{
    bool operator()(const pair<int,float> &a, const pair<int,float> &b)
    {
        return a.second < b.second;
    }
};


typedef priority_queue<pair<int,float>, vector<pair<int,float>>, DecreaseCmp > Decrease_queue;
typedef priority_queue<pair<int,float>, vector<pair<int,float>>, IncreaseCmp > Increase_queue;
typedef vector<pair<int,float>> Pair_vector;

bool myCmpfunc(pair<int,float> a, pair<int,float> b)
{
    return a.second > b.second;
}

struct _myclass {
    bool operator() (pair<int, Pair_vector> a, pair<int, Pair_vector> b)
    {
        return a.first < b.first;
    }
} myCmpClass;




template<typename T>
T removeExtra(T solution, int efK) {
    T result;
    while (!solution.empty())
    {
        result.push(solution.top());
        solution.pop();
        if (result.size() >= efK) break;
    }
    return result;
}

int deleteTop(vector<pair<int, float> >& vec, int size)
{
    // cout << vec.size() << " " << size << endl;
    vec.erase(vec.begin());
    pair<int, float> temp_pair;

    temp_pair.first = -1;
    temp_pair.second = numeric_limits<float>::min();
    vec.push_back(temp_pair);

    return 0;
}


int updateDecrease(vector<pair<int, float> >& vec, int size, pair<int, float> t)
{
    int i = size - 1;
    int j;
    if ((t.second < vec.back().second)) {
        return -1;
    }
    for (;;)
    {
        if (i == 0) break;
        j = i - 1;
        if (vec.at(j).second >= t.second) break;
        i = j;
    }
    j = size - 1;
    for (;;)
    {
        if (j == i) break;
        vec.at(j).second = vec.at(j-1).second;
        vec.at(j).first = vec.at(j-1).first;
        --j;
    }
    vec.at(i).second = t.second;
    vec.at(i).first = t.first;
    return 0;
}

// int updateIncrease(vector<pair<int, float> >& vec, int K_nei, pair<int, float> t)
// {
//     int i = K_nei - 1;
//     int j;
//     if ((t.second >= vec.back().second)) {
//         return -1;
//     }
//     for (;;)
//     {
//         if (i == 0) break;
//         j = i - 1;
//         if (vec.at(j).second < t.second) break;
//         i = j;
//     }
//     j = K_nei - 1;
//     for (;;)
//     {
//         if (j == i) break;
//         vec.at(j).second = vec.at(j-1).second;
//         vec.at(j).first = vec.at(j-1).first;
//         --j;
//     }
//     vec.at(i).second = t.second;
//     vec.at(i).first = t.first;
//     return 0;
// }

void searchGraph(FloatDataset& query, FloatDataset& data, int seed, int query_N, int N, int D, int K, int efK,
                vector<vector<int> >& neighbor, string out_search_file, string ground_path, float& tim, float& rec, const faiss::IndexIVFFlat* idx)
{
    int flag[N];
    vector<vector<int>> solu_vec;
    ofstream wfile(out_search_file);

    int nlist = idx->invlists->nlist;
    int n_bridges = 1000;   // !!!!!!!!!!!
    long* labels = new long [query_N * n_bridges];
    float * coarse_dis = new float [query_N * n_bridges];

    (dynamic_cast<faiss::MultiIndexQuantizer*>(idx->quantizer))->search_ip(query_N, query[0], n_bridges, coarse_dis, labels);  // ????????
    // add a function here : search_ip( XXXXX )
    // debug  there is still some problems with search_ip() function
    long total = query_N * n_bridges;
    for (int i = 0; i < total; ++i) {
        coarse_dis[i] = abs(coarse_dis[i]);
    }
    
    for (int j = 0; j < 5; ++j) {
        for (int i = 0; i < 3; ++i) {
            cout << "query_id "<< j << ":" << labels[j * n_bridges + i] << "," << coarse_dis[j * n_bridges + i] << "\t";
        }
        cout << endl;
    }
    

    auto start = std::chrono::high_resolution_clock::now();
    srand (seed);
    for (int i = 0; i < query_N; i++)
    {
        memset(flag, 0, sizeof(flag));
        vector<pair<int, float> > candidate(efK, make_pair(-1, numeric_limits<float>::min()));
        vector<pair<int, float> > solution(K, make_pair(-1, numeric_limits<float>::min()));

        int st_point = rand() % N;
        float dis = 0.0;
        for (int d = 0; d < D; d++)
        {
            dis += data[st_point][d] * query[i][d];
        }
        updateDecrease(candidate, efK, make_pair(st_point, dis));
        /*
        if (i < 10)
            cout << i << ", start point :  "<< dis << endl;
        */
        // bridge vector (-2)
        int bridge_offset = 0;
        updateDecrease(candidate, efK, make_pair(-2, coarse_dis[i * n_bridges + bridge_offset]));


        // int count_imi = 0;

        while (candidate.front().first != -1)
        {
            pair<int,float> top = candidate.front();
            deleteTop(candidate, efK);
            int curr = top.first;

            if (curr == -2) {
                int list_id = labels[i * n_bridges + bridge_offset];
                const long* ids = idx->invlists->get_ids(list_id);
                // count_imi += idx->invlists->list_size(list_id);
                for (int j = 0; j < idx->invlists->list_size(list_id); ++j) {
                    // if (flag[ids[j]] != 0) continue;
                    // flag[ids[j]] = 1;
                    float dis = 0.0;
                    for (int d = 0; d < D; d++)
                    {
                        dis += data[ids[j]][d] * query[i][d];
                    }

                    /*
                    if (i < 3 && count_imi < 200) {
                        cout << i << ", imi candidate dis: " << dis << endl;
                    }
                    */
                    updateDecrease(candidate, efK, make_pair(ids[j], dis));
                }
                bridge_offset++;
                if (bridge_offset < n_bridges) {
                    updateDecrease(candidate, efK, make_pair(-2, coarse_dis[i * n_bridges + bridge_offset]));
                }
            }
            else {

                if (flag[curr] != 0) continue;
                flag[curr] = 1;

                updateDecrease(solution, K, top);

                for (int j = 0; j < neighbor[curr].size(); j++)
                {
                    int nei = neighbor[curr][j];
                    if (flag[nei] == 0)
                    {
                        float dis = 0.0;
                        for (int d = 0; d < D; d++)
                        {
                            dis += data[nei][d] * query[i][d];
                        }
                        updateDecrease(candidate, efK, make_pair(nei, dis));
                    }
                }
            }
        }

        vector<int> sol;
        for (int ss = 0; ss < solution.size(); ss++)
        {
            wfile << solution[ss].first << " " << solution[ss].second << " ";
            sol.push_back(solution[ss].first);
        }
        wfile << endl;
        solu_vec.push_back(sol);
    }
    auto end = std::chrono::high_resolution_clock::now();
    cout << "efK: " << efK << " qt: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() / (double)query_N;
    //wResult << " " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() / (double)query_N;
    tim = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() / (double)query_N;
    wfile.close();

    //recall
    cout << "start recall" << endl;

    ifstream rgound(ground_path);
    int line, kk;
    rgound >> line >> kk;
    vector<vector<int> > ground_truth;
    for (int i = 0; i < query_N; i++)
    {
        int index, id;
        float dis;
        rgound >> index;
        vector<int> single_truth;
        for (int j = 0; j < kk; j++)
        {
            rgound >> id >> dis;
            single_truth.push_back(id);
        }
        ground_truth.push_back(single_truth);
    }
    rgound.close();

    //cout << ground_truth.size() << " " << ground_truth[0].size() << endl;
    //cout << K << endl;
    int count = 0;
    for (int i = 0; i < query_N; i++) {
        for (int j = 0; j < K; j++) {
            for (int y = 0; y < K; y++)
            {
                if (ground_truth[i][j] == solu_vec[i][y]) count++;
            }
        }
    }

    cout << " Recall: " << (1.0 * count / (query_N*K)) << endl;
    rec = (1.0 * count / (query_N*K));
}

Increase_queue find_neighbors(int i, int N, int D, int K_nei, FloatDataset& data)
{
    Increase_queue temp_que;
    for (int j = 0; j < N; j++)
    {
        if (j == i) continue;
        float inner_sum = 0.0;
        // Inner product
        for (int d=0; d<D; d++)
        {
            inner_sum += data[i][d]*data[j][d];
        }
            temp_que.push(make_pair(j,inner_sum));
            if (temp_que.size()>K_nei) temp_que.pop();
    }
    return temp_que;
}

int main(int argc, char *argv[])
{
   //netflix
    if (argc != 4) {
        cout << "Input error" << endl;
        return 0;
    }

    int N, D, query_N, K, efK;
    int K_nei;
    K_nei = atoi(argv[3]);
    cout << "K_nei: " << K_nei << endl;

    string data_path, query_path, ground_path, out_search_file, out_construct;

    string dataset(argv[1]);

    if (string(argv[1]) == "netflix")
    {
        N = 17770;
        D = 300;
        query_N = 1000;
        K = 10;
        efK = 64;
        data_path = "/data/liujie/gqr/data/netflix/netflix_base.fvecs";
        query_path = "/data/liujie/gqr/data/netflix/netflix_query.fvecs";
        ground_path = "/data/liujie/gqr/data/netflix/netflix_top10_product_groundtruth.lshbox";

        out_search_file = "./search_netflix" + to_string(efK) + ".txt";
        out_construct = "./construct_netflix_" + to_string(K_nei) + ".txt";
    }
    else if (string(argv[1]) == "yahoo")
    {
        N = 136736;
        D = 300;
        query_N = 1000;
        K = 10;
        efK = 32;
        data_path = "/data/liujie/gqr/data/yahoomusic/yahoomusic_base.fvecs";
        query_path = "/data/liujie/gqr/data/yahoomusic/yahoomusic_query.fvecs";
        ground_path = "/data/liujie/gqr/data/yahoomusic/yahoomusic_top10_product_groundtruth.lshbox";

        out_search_file = "./search_yahoo_" + to_string(efK) + ".txt";
        out_construct = "./construct_yahoo_" + to_string(K_nei) + ".txt";
    }
    else if (string(argv[1]) == "imagenet")
    {
        N = 2340373;
        D = 150;
        query_N = 1000;
        K = 10;
        efK = 64;
        data_path = "/data/liujie/gqr/data/imagenet/imagenet_base.fvecs";
        query_path = "/data/liujie/gqr/data/imagenet/imagenet_query.fvecs";
        ground_path = "/data/liujie/gqr/data/imagenet/imagenet_top10_product_groundtruth.lshbox";


        out_search_file = "./search_imagenet" + to_string(efK) + ".txt";
        out_construct = "./construct_imagenet_" + to_string(K_nei) + ".txt";
    }
    else return 0;


    FloatDataset data(D, N);
    data.loadFvecs(data_path);
    
    string faiss_outputfile = dataset + "_pq.index";
    //Construction KNN
    if (string(argv[2]) == "0")
    {
        cout << "start construction" << endl;
        
        //vector<pair<int, Increase_queue> > que;
        vector<pair<int, Pair_vector> > que;
        const clock_t begin_time = clock();

        #pragma omp parallel for
        for (int i = 0; i < N; i++)
        {
            vector<pair<int, float> > temp_que(K_nei, make_pair(-1, numeric_limits<float>::min()));
            
            for (int j = 0; j < N; j++)
            {
                if (j == i) continue;
                float inner_sum = 0.0;
                // Inner product
                for (int d = 0; d < D; d++)
                {
                    inner_sum += data[i][d] * data[j][d];
                }
                
                if (inner_sum >= temp_que.at(K_nei-1).second) {
                    updateDecrease(temp_que, K_nei, make_pair(j, inner_sum));
                }

            }
            #pragma omp critical
            {
                reverse(temp_que.begin(), temp_que.end());
                que.push_back(make_pair(i, temp_que));
            }

        }
        cout << "Time: " << float( clock () - begin_time ) / CLOCKS_PER_SEC;

        cout << "finish construction" << endl;
        sort(que.begin(), que.end(), myCmpClass);

        ofstream wfile(out_construct);
        for (int i = 0; i < N; i++)
        {
            wfile << que[i].first << " ";
            //while (!que[i].second.empty())
            for (int j = 0; j < que[i].second.size(); j++)
            {

                pair<int,float> t = que[i].second[j];
                wfile << t.first << " " << t.second << " ";
            }
            wfile << endl;
        }
        wfile.close();
    }


    //Search KNN
    if (string(argv[2]) == "1")
    {
        // load PQ indexes from index file
        size_t nhash = 2;
        size_t nbits_subq = int(log2 (N + 1) / 2);     // good choice in general; may be modified later
        size_t ncentroids = 1 << (nhash * nbits_subq);  // total # of centroids
        faiss::MultiIndexQuantizer quantizer(D, nhash, nbits_subq);

        // faiss::MetricType metric = faiss::METRIC_L2;
        faiss::MetricType metric = faiss::METRIC_L2;

        faiss::IndexIVFFlat index(&quantizer, D, ncentroids, metric);
        index.quantizer_trains_alone = true;
        // index.maintain_direct_map = true;

        index.train(N, data[0]);
        index.add(N, data[0]);

        // faiss::write_index(&index, faiss_outputfile.c_str());

        // auto idx = unique_ptr<faiss::IndexIVFFlat>(dynamic_cast<faiss::IndexIVFFlat*>(faiss::read_index(faiss_outputfile.c_str())));

        cout << "nlist : " << index.invlists->nlist << endl;

        vector<vector<int> > neighbor;

        cout << "Start searching" << endl;
        ifstream rfile(out_construct);
        if (!rfile.is_open())
        {
            cout << "no construction file" << endl;
            return 0;
        }
        for (int i = 0; i < N; i++)
        {
            vector<int> vec;
            int index;
            rfile >> index;
            if (index != i) cout << "Read error" << endl;
            for (int j = 0; j < K_nei; j++)
            {

                pair<int,float> t;
                rfile >> t.first >> t.second;
                vec.push_back(t.first);
            }
            neighbor.push_back(vec);
        }
        rfile.close();

        FloatDataset query(D, query_N);
        query.loadFvecs(query_path);

        int nlist = index.invlists->nlist;

        int tmp, cnt = 0;
        for (int i = 0; i < nlist; ++i) {
            // cout << labels[i] << "," << coarse_dis[i] << "\t";
            tmp = index.invlists->list_size(i);
            if (tmp > 0) {
                cnt++;
                // cout << i << ", list size : " << tmp << endl;
            }
        }
        cout << endl;

        cout << "non-empty count : " << cnt << endl;
        // cout << "list 0 size : " << idx->invlists->list_size(0) << endl;


        string result_file = "./qt_r_K_" + to_string(K_nei) + "_" + string(argv[1]) + ".txt";
        ofstream wResult(result_file);
        int seed = 10;
        //int st_point = rand() % N;
        wResult << "efK qt recall" << endl;
        vector<int> ser_list;
        //for (int ef = 1; ef < 100; ef+=5) ser_list.push_back(ef);
        //for (int ef = 100; ef < 300; ef+=20) ser_list.push_back(ef);
        //for (int ef = 300; ef < 800; ef+=50) ser_list.push_back(ef);
        //for (int ef = 800; ef < 2000; ef+=100) ser_list.push_back(ef);
        ser_list.push_back(800);
        // for (int ef = 12000; ef < 50000; ef+=5000) ser_list.push_back(ef);
        // for (int ef = 50000; ef < 100000; ef+=1000) ser_list.push_back(ef);

        float tim = 0.0;
        float rec = 0.0;

        for (int ef = 0; ef < ser_list.size(); ef++)
        {

            searchGraph(query, data, seed, query_N, N, D, K, ser_list[ef],
                        neighbor, out_search_file, ground_path, tim, rec, &index);
            wResult << ser_list[ef] << " " << tim << " " << rec << endl;
        }
        wResult.close();

    }
    cout << "finish searching" << endl;

    return 0;
}
