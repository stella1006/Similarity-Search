#include <iostream>
#include <vector>
#include <utility>
#include <cstdlib>
#include <unistd.h>
#include <time.h>
#include <math.h>
using namespace std;
#define N 1000
#define query_N 1000
#define D 2
#define K 5

float inner[N][N];
int flag[N][N];

bool cmp(pair<int,float> a, pair<int,float> b)
{
    return a.second > b.second;
}

int main()
{
    srand (time(NULL));
    vector<vector<float> > point;
    vector<vector<float> > scale_point;
    vector<vector<float> > qeury;

    for (int i = 0; i < N; i++)
    {
        vector<float> dimension;
        for (int d = 0; d < D; d++)
        {
            dimension.push_back(rand()%299 + 1);
            usleep(20);
        }
        usleep(30);
        int sum = 0;
        for (int d = 0; d < D; d++)
        {
            sum+=dimension[d]*dimension[d];
        }
        point.push_back(dimension);
        for (int d = 0; d < D; d++)
        {
            dimension[d]=dimension[d]*((rand()%10+90)*1.0/100)/sqrt(sum);
        }
        scale_point.push_back(dimension);
    }
    for (int i = 0; i < query_N; i++)
    {
        vector<float> dimension;
        for (int d = 0; d < D; d++)
        {
            dimension.push_back(rand()%299 + 1);
            usleep(20);
        }
        point.push_back(dimension);
    }


    //memset(inner, 0, sizeof(inner));
    //memset(flag, 0, sizeof(inner));

    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            for (int d=0; d<D; d++)
            {
                inner[i][j] += scale_point[i][d]*scale_point[j][d];
            }
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
            length[i] += scale_point[i][d]*scale_point[i][d];
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
        //if (InDegree[j]!=0)
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


    return 0;
}
