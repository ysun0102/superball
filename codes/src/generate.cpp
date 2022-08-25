//
// Created by Yao Sun on 2021/07/22.
// Last revised 2022/05/31.
//

#include "superball.hpp"



/////////////////////////////////////////////////////////////////


pthread_mutex_t data_mutex = PTHREAD_MUTEX_INITIALIZER;
#define DATA_LOCK() {do {pthread_mutex_lock(&data_mutex); } while (0);}
#define DATA_UNLOCK() {do {pthread_mutex_unlock(&data_mutex); } while (0);}




// compute an inequality, i.e. excludes points as many as possible.
flag_t Task::generate_inequality(Input *p_in) {

    assert(this->first.coeff.size() == 0);

    GRBEnv env = GRBEnv();

    if (p_in->n_workers != 1) {
        env.set(GRB_IntParam_LogToConsole, 0);
    }
    env.set(GRB_IntParam_Threads, p_in->n_threads_gurobi);
    env.set(GRB_IntParam_MIPFocus, 3);
    GRBModel model = GRBModel(env);


    // variables
    wind_t coeff_limit = wind_t((SC_COEFF - 100) / (this->reg.max_dist + 1));
    wind_t constant_limit = coeff_limit * (this->reg.max_dist + 1);

    const wind_t dim = p_in->dim;
    vector<GRBVar> z(dim + 1);
    for (wind_t j = 0; j < dim; j++) {
        z[j] = model.addVar(0, coeff_limit, 0, GRB_INTEGER);
    }
    z[dim] = model.addVar(1, constant_limit, 0, GRB_INTEGER);


#ifdef REGION_then_BORDER

    const wind_t weight_region = this->reg.border.size() + 1;
    const wind_t weight_border = 1;

#endif

#ifdef BORDER_then_REGION

    const wind_t weight_region = 1;
    const wind_t weight_border = this->reg.region.size() + 1;

#endif

#ifdef WEIGHTED_ONE

    const wind_t weight_region = 1;
    const wind_t weight_border = 1;

#endif



    // build model
    // Constraint 1: border
    GRBLinExpr obj = 0;
    vector<GRBVar> nb;
    for (auto diff : this->reg.border) {
        INEQ_L(diff, z, dim)
        model.addConstr(l >= 0);

        GRBVar b = model.addVar(0, 1, 0, GRB_BINARY);
        nb.push_back(b);
        l += -SC_COEFF * (1 - b);
        model.addConstr(l <= 0);
        obj += b * weight_border;
    }

    // Constraint 2: region points
    vector<GRBVar> rg;
    for (wind_t r_i = 0; r_i < this->reg.region.size(); r_i++) {
        p_t diff = this->reg.region[r_i];
        INEQ_L(diff, z, dim)

        GRBVar b = model.addVar(0, 1, 0, GRB_BINARY);
        rg.push_back(b);
        l += -SC_COEFF * (1 - b);
        model.addConstr(l <= -1);
        obj += b * weight_region;
    }

    // Constraint 3: exclude previous computed inequalities.
    for (auto ineq : this->ineq_set) {        
        if (ineq.noton_index.size()) {
            GRBLinExpr exclude = 0;           
            for (auto idx : ineq.noton_index) {
                exclude += nb[idx];
            }
            model.addConstr(exclude >= 1);
        }
    }
    for (auto ineq : this->ineq_set) {
        if (ineq.inc_index.size()) {
            GRBLinExpr exclude = 0;           
            for (auto idx : ineq.inc_index) {
                exclude += rg[idx];
            }
            model.addConstr(exclude >= 1);
        }
    }


    model.setObjective(obj, GRB_MAXIMIZE);

    double this_time = 0;
    for (wind_t durable = 0; durable < GENERATE_DURABLE; durable++) {
        model.set(GRB_DoubleParam_TimeLimit, GENERATE_TIMELIMIT);
        model.optimize();
        this_time += model.get(GRB_DoubleAttr_Runtime);
        if (model.get(GRB_IntAttr_SolCount) || model.get(GRB_IntAttr_Status) == GRB_INFEASIBLE) {
            break;
        }
    }
    this->time += this_time;

    if (model.get(GRB_IntAttr_SolCount) == 0) {
        return MODEL_INFESIBLE;
    }


    flag_t result = (model.get(GRB_IntAttr_Status) == GRB_TIME_LIMIT) ? MODEL_TIMELIMIT : MODEL_OPTIMAL;

    model.set(GRB_IntParam_SolutionNumber, 0);

    wind_t most_on = 0;
    for (auto b : nb) {
        if (round(b.get(GRB_DoubleAttr_Xn)) == 1) {
            most_on++;
        }
    }
    wind_t best_weight = 0;
    for (auto b : rg) {
        if (round(b.get(GRB_DoubleAttr_Xn))) {
            best_weight++;
        }
    }

    this->first.center = this->center;
    for (wind_t j = 0; j < dim + 1; j++) {
        this->first.coeff.push_back(round(z[j].get(GRB_DoubleAttr_Xn)));
    }
    this->first.update_inc_exc(this->reg);
    this->first.update_noton_on(this->reg);

    assert(this->first.exc_diff.size() >= best_weight);
    assert(this->first.on_diff.size() >= most_on);


    if (this->first.inc_index.size() == 0 && this->first.noton_index.size() == 0) {
        result = MODEL_INFESIBLE;
    }

    this->ineq_set.push_back(this->first);
    this->first.clear();

    return result;
}






void Task::superball_type_point(Input *p_in, Dyn_task *dyn) {

    // step 0: check prevous files
    char filename[255];
    sprintf(filename, "%s/%d.txt", p_in->ineqs_filename.c_str(), this->center);
    FILE *test = fopen(filename, "r");
    if (test) {
        PRINTF_STAMP("file '%s' exists, no needs of re-computation\n", filename);
        fclose(test);
        return;
    }

    // step 1: build diff_to_type table.
    this->reg.center = this->center;
    vector<flag_t> diff_to_type(1 << p_in->dim);
    this->reg.diff_to_type = diff_to_type;
    this->time = 0;

    // step 2: generate region.
    this->reg.generate(p_in);

    // step 3: compute GENERATE_NUM inequality
    for (wind_t i = 0; i < GENERATE_NUM; i++) {
        // compute an inequality, the inequality is stored in this->first.
        if (this->generate_inequality(p_in) == MODEL_INFESIBLE) {
            break;
        }
    }


    // step 6: output the inequalities in this->ineq_set to files
    FILE *fout = fopen(filename, "w");
    for (auto ineq : this->ineq_set) {
        ineq.write_detail(fout);
    }
    fclose(fout);

    
    cout << endl;
    PRINTF_STAMP("[%d]: center %d [%x] is DONE\n", \
        this->thread_idx, this->center, this->center);

    printf("\t\t#Region %3s %ld\n", ":", this->reg.region.size());
    printf("\t\t#Outter %3s %ld\n", ":", this->reg.outter.size());
    printf("\t\t#Border %3s %ld\n", ":", this->reg.border.size());
    printf("\t\t#mBorder %2s %ld\n", ":", this->reg.miniborder.size());
    printf("\t\tMax-dist %2s %d\n", ":", this->reg.max_dist);
    for (wind_t i = 0; i < this->ineq_set.size(); i++) {
        printf("\t\t%2d-th %5s exc %ld/%ld, on %ld/%ld\n", i, ":", \
            this->ineq_set[i].exc_diff.size(), this->reg.region.size(), \
            this->ineq_set[i].on_diff.size(), this->reg.border.size());
    }
    printf("\t\tTime %6s %5.3f\n\n", ":", this->time);

}







void *get_task(void *paraall) {

    parameters_t *par = (parameters_t *) paraall;
    Dyn_task *dyn = par->p_dyn;
    
    Task tk;
    bool i_am_working = true;

    while (true) {
        DATA_LOCK()
        if (dyn->task_array.size() == 0) {
            if (dyn->working_threads == 0) {               
                // all works are done
                DATA_UNLOCK()
                return NULL;

            } else {
                // still has workers
                if (i_am_working) {
                    dyn->working_threads--;
                    i_am_working = false;
                }
                DATA_UNLOCK()
                if (dyn->task_array.size()) {
                    sleep(SLEEP_TIME); // in seconds.    
                }
                continue;
            }
        }

        if (!i_am_working) {
            dyn->working_threads++;
            i_am_working = true;
        }

        // assign task
        tk = dyn->task_array[dyn->task_array.size() - 1];
        tk.task_id = dyn->task_array.size() - 1;
        dyn->task_array.pop_back();
        tk.thread_idx = par->thread_idx;
        DATA_UNLOCK()
        

        switch (tk.type) {
            case TASK_TYPE_POINT: {

                tk.superball_type_point(par->p_in, dyn);

                break;
            }
            default:
                break;
        }



        DATA_LOCK()
        dyn->n_done += 1;
        DATA_UNLOCK()
    }

    return NULL;

}








int main(int argc, const char* argv[]) { 


    Input in(argc, argv);

    PRINTF_STAMP("Start computing superball in parallel\n\n");

    Dyn_task dyn_task(in);


    // initialize workers
    pthread_t *threads = new pthread_t[in.n_workers];
    pthread_attr_t attr;
    pthread_attr_init(&attr);

    parameters_t *para = new parameters_t[in.n_workers];

    for (wind_t th = 0; th < in.n_workers; th++) {

        para[th].p_in = &in;
        para[th].p_dyn = &dyn_task;
        
        para[th].thread_idx = th;

        pthread_create(threads + th, &attr, get_task, (void *) (para + th));
    }

    for (wind_t th = 0; th < in.n_workers; th++) {
        pthread_join(threads[th], NULL);
    }

    delete threads;
    delete para;



    PRINTF_STAMP("All done ~~~~~\n\n");

    cout << getCurrentSystemTime() << endl;

}

