

#include <cstdlib>
#include <string>
using std::string;

void netflix_test(size_t, size_t, string, string, string, int M, int efConstruction);
int main(int argc, char** argv) {

    size_t vecsize = std::atoi(argv[1]);
    size_t vecdim = std::atoi(argv[2]);
    string lshboxPath = argv[3];
    string basePath = argv[4];
    string queryPath = argv[5];
    int M = 16;
    int efConstruction = 200;
    netflix_test(vecsize, vecdim, lshboxPath, basePath, queryPath, M, efConstruction);
    return 0;
};
