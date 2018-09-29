#include <iostream>
#include <fstream>
#include <queue>
#include <chrono>
#include <vector>
#include <utility>
#include "hnswlib/hnswlib.h"


#include <unordered_set>

using std::vector;
using std::pair;
using namespace std;
using namespace hnswlib;

/*
template <typename T>
void writeBinaryPOD(ostream& out, const T& podRef) {
    out.write((char*)&podRef, sizeof(T));
}

template <typename T>
static void readBinaryPOD(istream& in, T& podRef) {
    in.read((char*)&podRef, sizeof(T));
}*/
class StopW {
    std::chrono::steady_clock::time_point time_begin;
public:
    StopW() {
        time_begin = std::chrono::steady_clock::now();
    }

    float getElapsedTimeMicro() {
        std::chrono::steady_clock::time_point time_end = std::chrono::steady_clock::now();
        return (std::chrono::duration_cast<std::chrono::microseconds>(time_end - time_begin).count());
    }

    void reset() {
        time_begin = std::chrono::steady_clock::now();
    }

};

float test_approx(float *massQ, size_t vecsize, size_t qsize, HierarchicalNSW<float> &appr_alg, size_t vecdim,
                  vector<std::priority_queue<std::pair<float, labeltype >>> &answers, size_t k) {
    size_t correct = 0;
    size_t total = 0;
//#pragma omp parallel for
    for (int i = 0; i < qsize; i++) {

        std::priority_queue<std::pair<float, labeltype >> result = appr_alg.searchKnn(massQ + vecdim * i, k);
        std::priority_queue<std::pair<float, labeltype >> gt(answers[i]);
        unordered_set<labeltype> g;
        total += gt.size();
        while (gt.size()) {
            g.insert(gt.top().second);
            gt.pop();
        }
        while (result.size()) {
            if (g.find(result.top().second) != g.end())
                correct++;
            result.pop();
        }
    }
    return 1.0f * correct / total;
}

void test_vs_recall(float *massQ, size_t vecsize, size_t qsize, HierarchicalNSW<float> &appr_alg, size_t vecdim,
                    vector<std::priority_queue<std::pair<float, labeltype >>> &answers, size_t k) {
    //vector<size_t> efs = { 1,2,3,4,6,8,12,16,24,32,64,128,256,320 };//  = ; { 23 };
    vector<size_t> efs;
    for (int i = 10; i < 30; i+=10) {
        efs.push_back(i);
    }
    for (int i = 100; i < 1000; i += 200) {
        efs.push_back(i);
    }
    for (int i = 2000; i < 10000; i += 2000) {
        efs.push_back(i);
    }
    for (int i = 20000; i < 100000; i += 20000) {
        efs.push_back(i);
    }
    /*for (int i = 300; i <600; i += 20) {
        efs.push_back(i);
    }*/
    for (size_t ef : efs) {
        appr_alg.setEf(ef);
        StopW stopw = StopW();

        float recall = test_approx(massQ, vecsize, qsize, appr_alg, vecdim, answers, k);
        float time_us_per_query = stopw.getElapsedTimeMicro() / qsize;
        cout << ef << "\t" << recall << "\t" << time_us_per_query << " us\n";
        if (recall > 1.0) {
            cout << recall << "\t" << time_us_per_query << " us\n";
            break;
        }
    }
}

int loadFvecs(float*& data, int numItems, string inputPath) {
    ifstream fin(inputPath, ios::binary);
    if (!fin) {
        std::cout << "cannot open file " << inputPath << std::endl;
        assert(false);
    }
    int dimension;
    fin.read((char*)&dimension, 4);
    data = new float[numItems* dimension];
    fin.read((char*)data, sizeof(float) * dimension);

    int dim;
    for (int i = 1; i < numItems; ++i) {
        fin.read((char*)&dim, 4);
        assert(dim == dimension);
        fin.read((char*)(data + i * dimension), sizeof(float) * dimension);
    }
    fin.close();
    return dimension;
}

vector<std::priority_queue<std::pair<float, labeltype >>> loadLSHBOX(string inputPath) {
    vector<std::priority_queue<std::pair<float, labeltype >>> answers;

    ifstream fin(inputPath.c_str());
    unsigned numQueries;
    unsigned K;
    fin >> numQueries >> K;
    answers.resize(numQueries);

    unsigned qId;
    unsigned id;
    float dist;
    int index = 0;
    for (int q = 0; q < numQueries; ++q) { 
        fin >> qId;
        assert(qId == q);
        for (int i = 0; i < K; ++i) {
            fin >> id >> dist;
            answers[q].emplace(dist, id);
        }
    }
    fin.close();
    return answers;
}

void netflix_test(size_t vecsize, size_t vecdim, string lshboxPath, string basePath, string queryPath, int M = 16, int efConstruction = 200) {

    // L2Space l2space(vecdim);
    InnerProductSpace l2space(vecdim);

    vector<std::priority_queue<std::pair<float, labeltype >>> answers = loadLSHBOX(lshboxPath);
    size_t qsize = answers.size();
    size_t K = answers.front().size();


    // load base fvecs
    float *mass = NULL;
    size_t datadim = loadFvecs(mass, vecsize, basePath);
    assert(datadim == vecdim);

    // load query fvecs
    float *massQ = NULL;
    size_t querydim = loadFvecs(massQ, qsize, queryPath);
    assert(querydim == vecdim);


    HierarchicalNSW<float> appr_alg(&l2space, vecsize, M, efConstruction);

    cout << "Building index\n";
    StopW stopwb = StopW();
    for (int i = 0; i < 1; i++) {
        appr_alg.addPoint((void *) (mass + vecdim * i), (size_t) i);
    }
#pragma omp parallel for
    for (int i = 1; i < vecsize; i++) {
        appr_alg.addPoint((void *) (mass + vecdim * i), (size_t) i);
    }

    cout << "Index built, time=" << stopwb.getElapsedTimeMicro() * 1e-6 << "\n";
    string indexPath = "./";
    indexPath += std::to_string(vecsize) + "_" + std::to_string(vecdim) + "_"
        + std::to_string(M) + "_" + std::to_string(efConstruction) + "." + "hnswindex";
    appr_alg.saveIndex(indexPath);

    for (int i = 0; i < 1; i++)
        test_vs_recall(massQ, vecsize, qsize, appr_alg, vecdim, answers, K);
    cout << "finished";

    return;
}
