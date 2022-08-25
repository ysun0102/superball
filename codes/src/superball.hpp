//
// Created by Yao Sun on 2021/07/22.
// Last revised 2022/05/31.
//

#ifndef SUPERBALL_HPP
#define SUPERBALL_HPP



// ********************* inlcudes ************************

// standard lib
#include <limits.h>
#include <algorithm>
#include <iostream>
#include <assert.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <cstdio>
#include <bitset>
#include <vector>
#include <map>
#include <unordered_map>
#include <cmath>
#include <fstream>
#include <chrono>
#include <string.h>
#include <cstdlib>
#include "sys/time.h"

// pthread
#include <pthread.h>
#include <stdio.h>
#include <unistd.h>

// gubori
#include "gurobi_c++.h"

using namespace std;


// ********************** controllors **********************


//! controllors

// model basic
#define NUM_THREADS_WORKERS 32
#define NUM_THREADS_GUROBI 4

#define SLEEP_TIME 3

#define GENERATE_NUM 5
#define GENERATE_TIMELIMIT 30
#define GENERATE_DURABLE 5


// ********************** constants **********************


#define SC_COEFF 39999

// for region
#define REGION_INC 0
#define REGION_INC_BORDER 1
#define REGION_INC_MINIBORDER 2
#define REGION_EXC_CANNOT 3
#define REGION_EXC_INNER 4
#define REGION_EXC_OUTTER 5



#define MODEL_OPTIMAL 1
#define MODEL_INFESIBLE 2
#define MODEL_TIMELIMIT 3

#define TASK_TYPE_NONE 0
#define TASK_TYPE_POINT 1



#define FLAG_STRICT_INC 0
#define FLAG_ON 1
#define FLAG_EXC 2



// ! type redefinition
typedef int coeff_t;
typedef int wind_t;
typedef uint32_t p_t;
typedef vector<p_t> points_t;
typedef vector<wind_t> index_t;
typedef char flag_t;

#ifdef REGION_then_BORDER
#define ALG "REGION_then_BORDER"
#endif

#ifdef BORDER_then_REGION
#define ALG "BORDER_then_REGION"
#endif

#ifdef WEIGHTED_ONE
#define ALG "WEIGHTED_ONE"
#endif

#ifdef SIZE_then_STRENGTH
#define ALG "SIZE_then_STRENGTH"
#endif





//! useful functions
#define SB_MAX(a, b) ((a) > (b) ? (a) : (b))
#define SB_MIN(a, b) ((a) < (b) ? (a) : (b))


#define INEQ_L(diff, z, dim) \
    GRBLinExpr l = 0; \
    p_t tp = diff; \
    wind_t j = 0; \
    while (tp) { \
        if (tp & 1) { \
            l += z[j]; \
        } \
        tp >>= 1; \
        j++; \
    } \
    l -= z[dim];


/* convert timeval to miliseconds */
#define TIMEVAL2F(stamp) ((stamp).tv_sec * 1000.0 + (stamp).tv_usec / 1000.0)

/* get timestamp to the precision of miliseconds since the program starts */
inline double get_timestamp() {
    static double __init_stamp = -1;
    static struct timeval __cur_time;
    if (-1 == __init_stamp) {
        gettimeofday(&__cur_time, NULL);
        __init_stamp = TIMEVAL2F(__cur_time);
    }
    
    gettimeofday(&__cur_time, NULL);
    return ((TIMEVAL2F(__cur_time) - __init_stamp) / 1000.0);
}

/* print msg with timestamp */
#define PRINTF_STAMP(format, ...) \
do { \
    flockfile(stdout); \
    printf("%12.2f - ", get_timestamp()); \
    printf(format, ##__VA_ARGS__); \
    fflush(stdout); \
    funlockfile(stdout); \
} while(0)


inline string getCurrentSystemTime() {
    auto tt = chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    struct tm* ptm = localtime(&tt);
    char date[60] = { 0 };
    sprintf(date, "%d-%02d-%02d-%02d:%02d:%02d", (int)ptm->tm_year + 1900, (int)ptm->tm_mon + 1, (int)ptm->tm_mday, (int)ptm->tm_hour, (int)ptm->tm_min, (int)ptm->tm_sec);
    return string(date);
}

















class Input {

public:

    points_t inc;
    points_t exc;
    vector<flag_t> table; // 1 means inclusive; 0 means exclusive.

    wind_t dim;
    
    wind_t min_strength;

    wind_t n_workers;
    wind_t n_threads_gurobi;

    string points_filename;
    string ineqs_filename;
    
    vector<string> source_dirs;

    vector<vector<p_t>> weight_list;

    Input(int argc, const char * argv[]);
    void load_points();
    void calculate_weight_list();
    
}; 

class Region {

public:

    vector<flag_t> diff_to_type;
    index_t diff_to_index;

    p_t center;
    wind_t max_dist;

    points_t border;
    points_t miniborder;
    points_t region;
    points_t outter;

    void generate(Input *p_in);

};




class Inequality {
public:
    p_t center;

    vector<coeff_t> coeff;
    index_t inc_index;
    points_t on_diff;
    index_t noton_index;
    points_t exc_diff;

    flag_t point_notin_equ_normal(p_t diff);
    
    void update_inc_exc(Region &reg);
    void update_noton_on(Region &reg);
    
    void simplify();
    string tostring();
    wind_t zerocoeff();
    
    void display();
    void display_detail();
    void write(FILE *fout);
    void write_detail(FILE *fout);    
    bool read(FILE *fin, p_t center, wind_t dim);
    void clear();
};





class Dyn_task;

class Task {

public:

    wind_t type;
    p_t center;

    wind_t task_id;
    wind_t thread_idx;
    
    double time;

    Region reg;

    Inequality first;
    vector<Inequality> ineq_set;


    void superball_type_point(Input *p_in, Dyn_task *dyn);
    flag_t generate_inequality(Input *p_in);

};



class Dyn_task {

public:

    vector<Task> task_array;
    wind_t n_done;
    wind_t working_threads;

    Dyn_task(Input &in);

};


struct parameters_t {

    Input *p_in;
    Dyn_task *p_dyn;  
    wind_t thread_idx;

};




class Pool {

public:

    double time;

    wind_t cover_borders;

    vector<Inequality> candidate_set;
    
    index_t sel_index;
    vector<Inequality> selection;

    Pool(Input &in);
    
    void find_min_size(Input &in);
 
    void display();
    void write(Input &in);
    void clear();


};









// ********************** functions **********************



Input::Input(int argc, const char * argv[]) { 

    // points_filename, dim, n_workers, n_threads_gurobi
    if (argc < 3) {
        printf("Please input points_filename, dim, n_workers, n_threads_gurobi, and try agian.\n");
        exit(1);
    }

    this->points_filename = argv[1];
    this->dim = atoi(argv[2]);
    this->n_workers = (argc > 3) ? atoi(argv[3]) : NUM_THREADS_WORKERS;
    this->n_threads_gurobi = (argc > 4) ? atoi(argv[4]) : NUM_THREADS_GUROBI;
    

    {
        wind_t i = 0;
        while (this->points_filename[i] != '.') {
            this->ineqs_filename += this->points_filename[i++];
        }
    }
    this->ineqs_filename += "_";
    this->ineqs_filename += ALG;


    if (strcmp(ALG, "SIZE_then_STRENGTH") == 0) {
        if (argc < 7) {
            printf("Please input points_filename, dim, n_workers, n_threads_gurobi, min_strength, candidate_dir1, ...and try agian.\n");
            exit(1);
        }

        this->min_strength = atoi(argv[5]);
        for (wind_t i = 6; i < argc; i++) {
            this->source_dirs.push_back(argv[i]);
        }
        this->ineqs_filename += "_minstrength" + to_string(this->min_strength);

    } else {
        this->min_strength = -1;
        for (wind_t i = 5; i < argc; i++) {
            this->source_dirs.push_back(argv[i]);
        }
    }

    if (this->dim >= sizeof(p_t) * 8) {
        printf("Dimension is too large, %d vs. %ld.\n", this->dim, sizeof(p_t) * 8);
        exit(1);
    }



    printf("\n++++++++++++++++++++++++++++ parameter +++++++++++++++++++++++++++++++++\n\n");

    printf("Algorithm: %s\n", ALG);

    printf("Worker number: %d\n", this->n_workers);
    printf("Gurobi threads: %d\n\n", this->n_threads_gurobi);

    printf("Coefficient limits: %d\n\n", SC_COEFF);

    printf("GENERATE_NUM: %d\n", GENERATE_NUM);
    printf("GENERATE_TIMELIMIT: %d\n", GENERATE_TIMELIMIT);
    printf("GENERATE_DURABLE: %d\n\n", GENERATE_DURABLE);

    printf("Points filename: %s\n", this->points_filename.c_str());
    printf("Dimension: %d\n\n", this->dim);

    printf("Inequality dirname: %s\n\n", this->ineqs_filename.c_str());
    mkdir(this->ineqs_filename.c_str(), S_IRWXU);

    printf("Program: %s\n\n", argv[0]);

    if (this->min_strength > 0) {
        printf("Min_strength: %d\n\n", this->min_strength);
    }

    if (this->source_dirs.size()) {
        printf("Source dirs:\n");
        for (auto st : this->source_dirs) {
            printf("\t%s\n", st.c_str());
        }
        printf("\n\n");
    }

    cout << getCurrentSystemTime() << endl;

    printf("\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n\n");

    this->load_points();
    this->calculate_weight_list();
}





void Input::load_points() { 

    PRINTF_STAMP("loading points...\n");

    FILE *fin = fopen(("data/" + this->points_filename).c_str(), "r");
    this->inc.clear();
    p_t poi;
    while (fscanf(fin, "%x", &poi) > 0) {
        this->inc.push_back(poi);
    }
    fclose(fin);

    PRINTF_STAMP("to include points num: %ld\n", this->inc.size());


    // table
    vector<char> table(1 << this->dim);
    for (auto poi : this->inc) {
        table[poi] = 1;
    }
    this->table = table;


    // exclusive points    
    this->exc.clear();
    for (p_t i = 0; i < (1 << this->dim); i++) {
        if (table[i] == 0) {
            this->exc.push_back(i);
        }
    }
    PRINTF_STAMP("to exclude points num: %ld\n", this->exc.size());

    PRINTF_STAMP("points are loaded\n\n");

}



inline wind_t weight(p_t a) { 
    p_t c = a;
    wind_t res = 0;
    while (c != 0) {
        if (c & 1) {
            res += 1;
        }
        c >>= 1;
    }
    return res;
}




void Input::calculate_weight_list() {

    vector<uint64_t> table(1 << this->dim);
    for (p_t diff = 0; diff < (1 << this->dim); diff++) {
        table[diff] = uint64_t(diff) + (weight(diff) << this->dim);
    }
    sort(table.begin(), table.end());

    const p_t mask = (1 << this->dim) - 1;
    wind_t w = 0;
    points_t temp;
    for (auto poi : table) {
        if ((poi >> this->dim) == w) {
            temp.push_back(mask & poi);
        } else {
            this->weight_list.push_back(temp);
            temp.clear();
            w++;
            temp.push_back(mask & poi);
        }
    }
    assert(w == this->dim);
    this->weight_list.push_back(temp);

    PRINTF_STAMP("weight table is calculated\n\n");
    
}



// return true if not in.
flag_t Inequality::point_notin_equ_normal(p_t diff) {

    p_t tp = diff;
    wind_t j = 0;
    wind_t dim = this->coeff.size() - 1;
    coeff_t ans = -this->coeff[dim]; // in mormal form

    while (tp) {
        if (tp & 1) {
            ans += this->coeff[j];
        }
        tp >>= 1;
        j++;
    }

    if (ans > 0) {
        return FLAG_STRICT_INC;
    } else if (ans == 0) {
        return FLAG_ON;
    } else {
        return FLAG_EXC;
    }
}



void Inequality::update_inc_exc(Region &reg) { 
    assert(this->inc_index.size() == 0);
    assert(this->exc_diff.size() == 0);
    for (wind_t idx = 0; idx < reg.region.size(); idx++) {
        if (this->point_notin_equ_normal(reg.region[idx]) != FLAG_EXC) {
            this->inc_index.push_back(idx);
        } else {
            this->exc_diff.push_back(reg.region[idx]);
        }
    }
}


void Inequality::update_noton_on(Region &reg) {
    assert(this->noton_index.size() == 0);
    assert(this->on_diff.size() == 0);    
    for (wind_t idx = 0; idx < reg.border.size(); idx++) {
        if (this->point_notin_equ_normal(reg.border[idx]) != FLAG_ON) {
            this->noton_index.push_back(idx);
        } else {
            this->on_diff.push_back(reg.border[idx]);
        }
    }
}



void Inequality::simplify() {
    assert(this->coeff.size() > 1);

    coeff_t g = this->coeff[0];

    for (wind_t j = 1; j < this->coeff.size(); j++) {
        g = __gcd(g, this->coeff[j]);
    }
    g = g > 0 ? g : -g;
    assert(g > 0);
    if (g == 1) {
        return;
    } else {
        for (wind_t j = 0; j < this->coeff.size(); j++) {
            this->coeff[j] /= g;
        }
    }
}




string Inequality::tostring() {
    string result;
    coeff_t cnst = 0;
    for (wind_t j = 0; j < this->coeff.size() - 1; j++) {
        if (j) {
            result += "+";
        }
        if ((this->center >> j) & 1) {
            result += to_string(-this->coeff[j]);
            cnst += this->coeff[j];
        } else {
            result += to_string(this->coeff[j]);
        }
    }
    result += ">=" + to_string(cnst - this->coeff[this->coeff.size() - 1]);
    
    return result;
}


wind_t Inequality::zerocoeff() {
    wind_t res = 0;
    for (wind_t j = 0; j < this->coeff.size() - 1; j++) {
        if (this->coeff[j] == 0) {
            res++;
        }
    }   
    return res;
}





void Inequality::display() {
    coeff_t cnst = 0;
    for (wind_t j = 0; j < this->coeff.size() - 1; j++) {
        if ((this->center >> j) & 1) {
            printf("%5d, ", -this->coeff[j]);
            cnst += this->coeff[j];
        } else {
            printf("%5d, ", this->coeff[j]);
        }
    }
    printf("%5d [exc %5ld, on %5ld]", cnst - this->coeff[this->coeff.size() - 1], this->exc_diff.size(), this->on_diff.size());
    cout << endl;
}




void Inequality::display_detail() {
    coeff_t cnst = 0;
    for (wind_t j = 0; j < this->coeff.size() - 1; j++) {
        if ((this->center >> j) & 1) {
            printf("%d, ", -this->coeff[j]);
            cnst += this->coeff[j];
        } else {
            printf("%d, ", this->coeff[j]);
        }
    }
    printf("%d [on %ld:", cnst - this->coeff[this->coeff.size() - 1], this->on_diff.size());
    for (auto diff : this->on_diff) {
        printf(" %x", diff);
    }
    printf("], {exc %ld:", this->exc_diff.size());
    for (auto diff : this->exc_diff) {
        printf(" %x", diff);
    }
    printf("}");
    cout << endl;
}

void Inequality::write(FILE *fout) {
    coeff_t cnst = 0;
    for (wind_t j = 0; j < this->coeff.size() - 1; j++) {
        if ((this->center >> j) & 1) {
            fprintf(fout, "%d, ", -this->coeff[j]);
            cnst += this->coeff[j];
        } else {
            fprintf(fout, "%d, ", this->coeff[j]);
        }
    }
    fprintf(fout, "%d\n", cnst - this->coeff[this->coeff.size() - 1]);
}

void Inequality::write_detail(FILE *fout) {
    coeff_t cnst = 0;
    for (wind_t j = 0; j < this->coeff.size() - 1; j++) {
        if ((this->center >> j) & 1) {
            fprintf(fout, "%d, ", -this->coeff[j]);
            cnst += this->coeff[j];
        } else {
            fprintf(fout, "%d, ", this->coeff[j]);
        }
    }
    fprintf(fout, "%d::[%ld:", cnst - this->coeff[this->coeff.size() - 1], this->on_diff.size());
    for (auto diff : this->on_diff) {
        fprintf(fout, " %x", diff);
    }
    fprintf(fout, "]:::{%ld:", this->exc_diff.size());
    for (auto diff : this->exc_diff) {
        fprintf(fout, " %x", diff);
    }
    fprintf(fout, "}\n");
}




bool Inequality::read(FILE *fin, p_t center, wind_t dim) {

    assert(this->coeff.size() == 0);
    assert(this->exc_diff.size() == 0);

    this->center = center;

    coeff_t cnst = 0;
    coeff_t c;
    if (fscanf(fin, "%d, ", &c) > 0) {
        if (c < 0) {
            assert(center & p_t(1));
            c = -c;
            cnst += c;
        }
        this->coeff.push_back(c);
    } else {
        return false;
    }

    for (wind_t j = 1; j < dim; j++) {
        if (fscanf(fin, "%d, ", &c) <= 0) {
            exit(1);
        }
        if (c < 0) {
            assert(center & (p_t(1) << j));
            c = -c;
            cnst += c;
        }

        this->coeff.push_back(c);
    }
    wind_t len;
    if (fscanf(fin, "%d::[%d:", &c, &len) <= 0) {
        exit(1);
    }
    c = -c + cnst;
    assert(c > 0);
    this->coeff.push_back(c);

    p_t on;
    for (wind_t i = 0; i < len; i++) {
        if (fscanf(fin, " %x", &on) <= 0) {
            exit(1);
        }
        this->on_diff.push_back(on);
    }
    
    if (fscanf(fin, "]:::{%d:", &len) <= 0) {
        exit(1);
    }

    p_t ex;
    for (wind_t i = 0; i < len - 1; i++) {
        if (fscanf(fin, " %x", &ex) <= 0) {
            exit(1);
        }
        this->exc_diff.push_back(ex);
    }
    if (fscanf(fin, " %x}\n", &ex) <= 0) {
        exit(1);
    }
    this->exc_diff.push_back(ex);

    return true;
}



void Inequality::clear() {
    this->coeff.clear();
    
    this->noton_index.clear();
    this->on_diff.clear();

    this->inc_index.clear();
    this->exc_diff.clear();
}




Dyn_task::Dyn_task(Input &in) {

    for (wind_t idx = 0; idx < in.exc.size(); idx++) {

        Task cur_task;

        cur_task.type = TASK_TYPE_POINT;
        cur_task.center = in.exc[idx];

        cur_task.task_id = idx;
        cur_task.thread_idx = -1;

        this->task_array.push_back(cur_task);
    }

    this->n_done = 0;
    this->working_threads = in.n_workers;

    PRINTF_STAMP("initialized %ld tasks, type: TASK_TYPE_POINT\n\n", this->task_array.size());
}





void Region::generate(Input *p_in) {

    for (auto poi : p_in->exc) {
        this->diff_to_type[this->center ^ poi] = REGION_EXC_CANNOT;
    }

    // deal with weight, can break is no candidates.
    this->max_dist = 0;
    for (wind_t w = 0; w <= p_in->dim; w++) {
        bool has_outter = false;

        for (p_t diff : p_in->weight_list[w]) {
            points_t subsetpoints;
            wind_t j = 0;
            p_t td = diff;
            while (td) {
                if (td & p_t(1)) {
                    subsetpoints.push_back(diff ^ (p_t(1) << j));
                }
                td >>= 1;
                j++;
            }

            if (this->diff_to_type[diff] == REGION_INC) {
                this->diff_to_type[diff] = REGION_INC_MINIBORDER;
                for (auto subspoi : subsetpoints) {
                    if (this->diff_to_type[subspoi] <= REGION_EXC_CANNOT) {
                        this->diff_to_type[diff] = REGION_INC;
                        break;
                    }
                }

                if (this->diff_to_type[diff] == REGION_INC) {                    
                    for (auto subspoi : subsetpoints) {
                        if (this->diff_to_type[subspoi] > REGION_EXC_CANNOT) {
                            this->diff_to_type[diff] = REGION_INC_BORDER;
                            break;
                        }
                    }
                }
            } else {
                assert(this->diff_to_type[diff] == REGION_EXC_CANNOT);
                this->diff_to_type[diff] = REGION_EXC_OUTTER;
                for (auto subspoi : subsetpoints) {
                    if (this->diff_to_type[subspoi] <= REGION_EXC_CANNOT) {
                        this->diff_to_type[diff] = REGION_EXC_CANNOT;
                        break;
                    }
                }

                if (this->diff_to_type[diff] == REGION_EXC_OUTTER) {
                    has_outter = true;
                    for (auto subspoi : subsetpoints) {
                        assert(this->diff_to_type[subspoi] >= REGION_EXC_INNER);
                        this->diff_to_type[subspoi] = REGION_EXC_INNER;
                    }
                }
            }
        }

        if (!has_outter) {
            this->max_dist = w - 1;
            assert(this->max_dist >= 0);
            break;
        }
    }


    // update border and region
    this->miniborder.clear();
    this->region.clear();
    this->outter.clear();
    for (wind_t w = this->max_dist + 1; w >= 0; w--) {
        for (p_t diff : p_in->weight_list[w]) {
            switch (this->diff_to_type[diff]) {
                case REGION_INC_BORDER: {
                    this->border.push_back(diff);
                    break;
                }
                case REGION_INC_MINIBORDER: {
                    this->border.push_back(diff);
                    this->miniborder.push_back(diff);
                    break;
                }
                case REGION_EXC_INNER: {
                    this->region.push_back(diff);
                    break;
                }
                case REGION_EXC_OUTTER: {
                    this->region.push_back(diff);
                    this->outter.push_back(diff);
                    break;
                }
                default: {                    
                    break;
                }
            }
        }
    }
}






Pool::Pool(Input &in) {

    this->candidate_set.clear();
    char filename[255];
    assert(in.dim > 0);

    wind_t missing = 0;
    wind_t count = 0;
    wind_t ineq_count = 0;

    unordered_map<string, Inequality> umap;

    for (auto dir_name : in.source_dirs) {
        for (auto center : in.exc) {
            sprintf(filename, "%s/%d.txt", dir_name.c_str(), center);
            FILE *fin = fopen(filename, "r");
            if (!fin) {
                missing++;
                count++;
                continue;
            }

            if (cout && count % 10000 == 0) {
                PRINTF_STAMP("read %d / %ld files...\n", count, in.exc.size());
            }
            count++;

            while (true) {
                Inequality ineq;
                if (!ineq.read(fin, center, in.dim)) {
                    break;
                }
                ineq.simplify();
                umap[ineq.tostring()] = ineq;
                ineq_count++;
            }
            fclose(fin);
        }
    }


    for (auto it : umap) {
        this->candidate_set.push_back(it.second);
    }

    
    this->cover_borders = 0;

    PRINTF_STAMP("Loaded %d candidate inequalities from %ld centers, simplified to %ld, the other %d are missing\n\n", ineq_count, in.exc.size() - missing, this->candidate_set.size(), missing);
}






void Pool::display() {

    printf("#inequalities: %ld, cover borders: %d\n", this->selection.size(), this->cover_borders);
    printf("Selection: \n");
    for (wind_t i = 0; i < this->sel_index.size(); i++) {
        printf("eq%4d: ", this->sel_index[i]);
        this->selection[i].display();
    }
    cout << endl;
}


void Pool::write(Input &in) {
    char filename[255];
    sprintf(filename, "%s/selection_%ld.txt", in.ineqs_filename.c_str(), this->sel_index.size());

    FILE *fout = fopen(filename, "w");

    fprintf(fout, "%ld", this->sel_index.size());
    for (auto idx : this->sel_index) {
        fprintf(fout, " %d", idx);
    }
    fprintf(fout, "\n");

    for (auto ineq : this->selection) {
        ineq.write(fout);
    }
    fclose(fout);
}



void Pool::clear() {

    this->cover_borders = 0;

    this->sel_index.clear();
    this->selection.clear();

}











#endif /* SUPERBALL_HPP */
