#include <iostream>
#include <data.h>
#include <vector>
#include <utility>
#include <queue>
#include <algorithm>
using namespace std;

//imagenet
#define D 150
#define N 2340373
//netflix
// #define N 17770
// #define D 300
#define query_N 1000
#define K 10

//float InnerProduct[N][N];
vector<pair<int,float> > Norm_length;
long long ranking[N];
long long solution[1000];
// long long query_N;

typedef Dataset<float> FloatDataset;

bool cmp(pair<int,float> a, pair<int,float> b)
{
    return (b.second-a.second < 1e-6);
}

struct cmp_greater{
    bool operator() ( pair<int,float> a, pair<int,float> b){
        return (b.second-a.second < 1e-6);
    }
};

vector<std::priority_queue<std::pair<float, int >>> loadLSHBOX(string inputPath) {
    vector<std::priority_queue<std::pair<float, int >>> answers;

    ifstream fin(inputPath.c_str());
    unsigned numQueries;
    unsigned KK;
    fin >> numQueries >> KK;
    answers.resize(numQueries);

    unsigned qId;
    unsigned id;
    float dist;
    int index = 0;
    for (int q = 0; q < numQueries; ++q) {
        fin >> qId;
        assert(qId == q);
        for (int i = 0; i < KK; ++i) {
            fin >> id >> dist;
            answers[q].emplace(dist, id);
        }
    }
    fin.close();
    return answers;
}

int main()
{
    FloatDataset data(D, N);
    //FloatDataset query(D, query_N);
    vector<std::priority_queue<std::pair<float, int >>> query_ls;

    //imagenet
    query_ls = loadLSHBOX("/data/jinfeng/project/github/gqr/data/imagenet/imagenet_product_groundtruth.lshbox");
    data.loadFvecs("/data/jinfeng/project/github/gqr/data/imagenet/imagenet_base.fvecs");
    //query.loadFvecs("/data/jinfeng/project/github/gqr/data/imagenet/imagenet_query.fvecs");

    //netflix
    //query_ls = loadLSHBOX("/data/jinfeng/project/github/gqr/data/netflix/netflix_product_groundtruth.lshbox");
    //data.loadFvecs("/data/jinfeng/project/github/gqr/data/netflix/netflix_base.fvecs");
    //query.loadFvecs("/data/jinfeng/project/github/gqr/data/netflix/netflix_query.fvecs");

    cout << "Calaulating the norm......" << endl;
    // for (int i = 0; i < 10; i++)
    // {
    //     cout << data[i][0] << " " << data[i][1] << endl;
    // }
    for (long long i = 0; i < N; i++)
    {
        float norm = 0.0;
        for (int k = 0; k < D; k++)
        {
            norm += data[i][k]*data[i][k];
        }
        Norm_length.push_back(make_pair(i, norm));
    }
    sort(Norm_length.begin(), Norm_length.end(), cmp);

    cout << "Checking sort norm..." << endl;
    for (int i = 0; i < 10 ; i++)
    {
        cout << Norm_length[i].first << " " << Norm_length[i].second << endl;
    }
    cout << endl;

    cout << "Ranking......" << endl;
    #pragma omp parallel for
    for (long long i = 0; i < N; i++)
    {
        ranking[Norm_length[i].first] = i;
    }

    cout << "Finding solution to query......" << endl;
    cout << query_ls.size() << endl;
    for (int i = 0; i < query_N; i++)
    {
        //float record_min = 99999.0;

        // priority_queue<pair<int, float>, vector<pair<int, float> >, cmp_greater> min_que;
        // for (long long j = 0; j < N; j++)
        // {
        //     float inner = 0.0;
        //     for (int k = 0; k < D; k++)
        //     {
        //         inner += query[i][k]*data[j][k];
        //     }
        //     min_que.push(make_pair(j,inner));
        //     if (min_que.size()>K)
        //     {
        //         min_que.pop();
        //     }
        // }
        //get the last Element that rank topest
        long long Lowest_rank = 0;
        while (query_ls[i].size()>0) {
            //cout << min_que.top().second << " ";
            if (ranking[query_ls[i].top().second] > Lowest_rank)
            {
                Lowest_rank = ranking[query_ls[i].top().second];
            }
            query_ls[i].pop();
        }
        //cout << endl;
        solution[i] = Lowest_rank;
    }
    cout << query_N << endl;
    for (int i = 0; i < query_N; i++)
    {
        //cout << query_N << " ";
        //long long temp = solution[i];
        cout << i << ": " <<  solution[i] << " " << solution[i]*100.0/N << '%' <<  endl;
        //solution[i] << " " << solution[i]*100.0/N << '%' <<
        //cout << " " << data[temp][147] << " " << data[temp][148] << " " << data[temp][149];
        //cout << " " << query[i][147] << " " << query[i][148] << " " << query[i][149] << endl;
    }

    ofstream wfile("./output.csv");
    if (wfile.is_open())
    {
        for (int i = 0; i < query_N; i++)
        {
            //int temp = solution[i];
            wfile << solution[i] << ',' << solution[i]*100.0/N << endl;
        }
    }
    wfile.close();


    return 0;
}
