#include <iostream>
#include <data.h>
#include <vector>
#include <utility>
#include <queue>
#include <algorithm>
using namespace std;

//netflix
#define N 17770
#define D 300
#define query_N 1000
#define K 10

vector<vector<float> > inner(N);
vector<vector<int> > flag(N);

typedef Dataset<float> FloatDataset;

bool cmp(pair<int,float> a, pair<int,float> b)
{
    return a.second > b.second;
}

int main()
{
    FloatDataset data(D, N);
    data.loadFvecs("/data/jinfeng/project/github/gqr/data/netflix/netflix_base.fvecs");

    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            float inner_sum = 0.0;
            for (int d=0; d<D; d++)
            {
                inner_sum += data[i][d]*data[j][d];
            }
            inner[i].push_back(inner_sum);
        }
    }

    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            flag[i].push_back(0);
        }
    }

    for (int i = 0; i < N; i++)
    {
        vector<pair<int,float> > max_;
        for (int j = 0; j < N; j++)
        {
            max_.push_back(make_pair(j,inner[i][j]));
        }
        sort(max_.begin(), max_.end(), cmp);
        for (int k = 0; k < K; k++)
        {
            flag[i][max_[k].first] = 1;
        }
    }

    int InDegree[N];
    float length[N];
    for (int i = 0; i < N; i++)
    {
        length[i]=0;
        for (int d=0; d<D; d++)
            length[i] += data[i][d]*data[i][d];
    }

    for (int j = 0; j < N; j++)
    {
        int count=0;
        for (int i = 0; i < N; i++)
        {
            if (flag[i][j]==1) count++;
        }
        InDegree[j] = count;
    }

    for (int j = 0; j < N; j++)
    {
        if (InDegree[j]!=0)
            cout << j << " length: " << length[j] << " indegree: "<< InDegree[j] << endl;
            //<< " Point: "<< point[j][0] << " "<<  point[j][1]<< endl;;
    }


    cout << endl;
    int hasInDegree=0;
    for (int i = 0; i < N; i++)
    {
        if (InDegree[i] > 0) hasInDegree++;
    }

    cout << N << ": " << hasInDegree << endl;

    ofstream wfile("./eKNN_indegree.csv");
    if (wfile.is_open())
    {
        wfile << N << ',' << hasInDegree << endl;
        for (int j = 0; j < N; j++)
        {
            if (InDegree[j]!=0)
                wfile << j << ',' << length[j] << ',' << InDegree[j] << endl;
                //<< " Point: "<< point[j][0] << " "<<  point[j][1]<< endl;;
        }
    }
    wfile.close();

    return 0;
}
