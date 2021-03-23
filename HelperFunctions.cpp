

#include <unordered_map>
#include <unordered_set>
#include "HelperFunctions.h"


void ltrim(std::string &_str, const std::string &_chars = "\t\n\v\f\r ") {
    _str.erase(0, _str.find_first_not_of(_chars));
}

void rtrim(std::string &_str, const std::string &_chars = "\t\n\v\f\r ") {
    _str.erase(_str.find_last_not_of(_chars) + 1);
}

void trim(std::string &_str, const std::string &_chars = "\t\n\v\f\r ") {
    ltrim(_str, _chars);
    rtrim(_str, _chars);
}

/*
 *  Get the index of next unblank char from a string.
 */
unsigned int getNextChar(const char *str) {
    unsigned int rtn = 0;
    // Jump over all blanks
    for (; str[rtn] == ' '; ++rtn);
    return rtn;
}

/*
 *  Get next word from a string.
 */
std::string getNextWord(const char *str) {
    // Jump over all blanks
    std::string rtn(str);
    trim(rtn);
    return rtn;
}

double getCurrentTime() {
    long long time = std::chrono::duration_cast<std::chrono::microseconds>(
            std::chrono::steady_clock::now().time_since_epoch()).count();
    return static_cast<double>(time) / 1000000.0;
//    return clock() / (double) CLOCKS_PER_SEC;
}


Param param;

Param parseArgs(int nargs, char **args) {
    std::vector<std::string> Algorithms{
            "PowItr", "PowForPush", "SpeedPPR", "FwdPush", "GROUND_TRUTH", "BUILD_INDEX",
            "CLEAN_GRAPH", "GEN_QUERY"
    };

    Param rtn;
    std::unordered_set<std::string> algo_set(Algorithms.begin(), Algorithms.end());

    printf("%s\n", std::string(110, '=').c_str());
    printf("The input parameters are:\n");
    for (unsigned int cnt = 1; cnt < nargs;) {
        char *arg = args[cnt++];
        if (cnt == nargs) {
            printf("Unknown Parameters.\n");
            exit(0);
        }
        unsigned int i = getNextChar(arg);
        if (arg[i] != '-') {
            printf("Unknown Parameters.\n");
            exit(0);
        }
        std::string para = getNextWord(arg + i + 1);
        printf("-%s\n", para.c_str());
        arg = args[cnt++];
        if (para == "algo") {
            rtn.algorithm = std::string(arg);
            std::cout << "Algorithm Parameter: " << rtn.algorithm << std::endl;
            if (algo_set.find(rtn.algorithm) == algo_set.end()) {
                printf("Algorithm Unknown.\n");
                exit(0);
            }
        } else if (para == "graph") {
            rtn.graph_file = getNextWord(arg);
            printf("Input Graph File: %s\n", rtn.graph_file.c_str());
        } else if (para == "graph_binary") {
            rtn.graph_binary_file = getNextWord(arg);
            printf("Input Graph Binary File: %s\n", rtn.graph_binary_file.c_str());
        } else if (para == "with_idx") {
            auto option = getNextWord(arg);
            if (option == "yes") {
                rtn.with_idx = true;
            } else if (option == "no") {
                rtn.with_idx = false;
            } else {
                printf("Unknown option -%s!\n", option.c_str());
                exit(0);
            }
            std::cout << "With Index: " << rtn.with_idx << "\n";
        } else if (para == "is_undirected") {
            auto option = getNextWord(arg);
            if (option == "yes") {
                rtn.is_undirected_graph = true;
            } else if (option == "no") {
                rtn.is_undirected_graph = false;
            } else {
                printf("Unknown option -%s!\n", option.c_str());
                exit(0);
            }
            std::cout << "Is Undirected: " << rtn.is_undirected_graph << "\n";
        } else if (para == "meta") {
            rtn.meta_file = getNextWord(arg);
            printf("Input Meta File: %s\n", rtn.meta_file.c_str());
        } else if (para == "query") {
            rtn.query_file = getNextWord(arg);
            printf("Input Query File: %s\n", rtn.query_file.c_str());
        } else if (para == "answer_folder") {
            rtn.answer_folder = getNextWord(arg);
            printf("Answer Folder: %s\n", rtn.answer_folder.c_str());
        } else if (para == "output_folder") {
            rtn.output_folder = getNextWord(arg);
            printf("Output Folder: %s\n", rtn.output_folder.c_str());
        } else if (para == "estimation_folder") {
            rtn.output_estimations = true;
            rtn.estimation_folder = getNextWord(arg);
            printf("Estimation Folder: %s\n", rtn.estimation_folder.c_str());
        } else if (para == "index_file") {
            rtn.index_file = getNextWord(arg);
            printf("Output Index File: %s\n", rtn.index_file.c_str());
        } else if (para == "epsilon") {
            rtn.epsilon = std::stod(getNextWord(arg));
            printf("Epsilon: %.12f\n", rtn.epsilon);
        } else if (para == "alpha") {
            rtn.alpha = std::stod(getNextWord(arg));
            printf("Alpha: %.12f\n", rtn.alpha);
        } else if (para == "l1_error") {
            rtn.l1_error = std::stod(getNextWord(arg));
            rtn.specified_l1_error = true;
            printf("L1 ERROR: %.12f\n", rtn.l1_error);
        } else if (para == "query_size") {
            rtn.query_size = std::stoul(getNextWord(arg));
            printf("#Queries: %d\n", rtn.query_size);
        }  else {
            printf("Unknown option -%s!\n\n", para.c_str());
            exit(0);
        }
    }
    printf("%s\n", std::string(110, '-').c_str());
    /*****************************************************************************/
    return rtn;
}