//
// Created by Yao Sun on 2021/07/22.
// Last revised 2022/05/31.
//

#include "superball.hpp"



/////////////////////////////////////////////////////////////////



void Pool::find_min_size(Input &in) {

    GRBEnv env = GRBEnv();
    env.set(GRB_IntParam_Threads, in.n_threads_gurobi);
    env.set(GRB_IntParam_MIPFocus, 3);
    GRBModel model = GRBModel(env);


    PRINTF_STAMP("Start building model for %ld candidates\n\n", this->candidate_set.size());


    // all contained points
    vector<GRBVar> ineq_var;
    vector<GRBVar> on_var;

    GRBLinExpr obj = 0;
    for (auto ineq : this->candidate_set) {        
        GRBVar v = model.addVar(0, 1, 0, GRB_BINARY);
        ineq_var.push_back(v);
        obj += v * (in.inc.size() + 1);
    }

    vector<GRBLinExpr> exc_point_expr(1 << in.dim);
    vector<GRBLinExpr> on_point_expr(1 << in.dim);
    vector<bool> on_point_flag(1 << in.dim);
    
    for (auto ic : in.inc) {
        GRBVar v = model.addVar(0, 1, 0, GRB_BINARY);
        on_var.push_back(v);
        obj -= v;
    }

    for (wind_t inq_i = 0; inq_i < this->candidate_set.size(); inq_i++) {
        if (inq_i && inq_i % 10000 == 0) {
            PRINTF_STAMP("%d inequalities are done\n", inq_i);
        }
        for (auto diff : this->candidate_set[inq_i].exc_diff) {
            exc_point_expr[diff ^ this->candidate_set[inq_i].center] += ineq_var[inq_i];
        }
        for (auto diff : this->candidate_set[inq_i].on_diff) {
            on_point_expr[diff ^ this->candidate_set[inq_i].center] += ineq_var[inq_i];
            on_point_flag[diff ^ this->candidate_set[inq_i].center] = true;
        }
    }

    for (auto ex : in.exc) {
        model.addConstr(exc_point_expr[ex] >= 1);
    }

    for (wind_t i = 0; i < in.inc.size(); i++) {
        if (on_point_flag[in.inc[i]]) {
            model.addConstr(on_point_expr[in.inc[i]] >= in.min_strength);
            model.addConstr(on_point_expr[in.inc[i]] >= on_var[i]);
        } else {
            model.addConstr(on_var[i] == 0);
        }
    }


    if (this->sel_index.size() > 0) {
        model.addConstr(obj <= this->sel_index.size() - 1);
    }


    PRINTF_STAMP("Building is done\n\n");


    // start optimization
    model.setObjective(obj, GRB_MINIMIZE);

    model.optimize();
    this->time += model.get(GRB_DoubleAttr_Runtime);

    if (model.get(GRB_IntAttr_Status) == GRB_INFEASIBLE) {
        return;
    }
    
    assert(model.get(GRB_IntAttr_SolCount) > 0);
    model.set(GRB_IntParam_SolutionNumber, 0);

    this->clear();
    on_point_flag.clear();
    for (wind_t i = 0; i < ineq_var.size(); i++) {
        if (round(ineq_var[i].get(GRB_DoubleAttr_Xn)) == 1) {
            this->sel_index.push_back(i);
            this->selection.push_back(this->candidate_set[i]);

            for (auto diff : this->candidate_set[i].on_diff) {
                on_point_flag[diff ^ this->candidate_set[i].center] = true;
            }
        }
    }
    for (auto ic : in.inc) {
        if (on_point_flag[ic]) {
            this->cover_borders++;
        }
    }

    this->write(in);

    PRINTF_STAMP("best result %ld, cover borders %d\n\n", this->sel_index.size(), this->cover_borders);

}


int main(int argc, const char* argv[]) {

    Input in(argc, argv);

    PRINTF_STAMP("Select inequalities from candidate set\n\n");

    // load candidate inequalities.
    Pool pool(in);

    // find the minimial number of generators.
    pool.find_min_size(in);


    PRINTF_STAMP("All done ~~~~~\n\n");

    cout << getCurrentSystemTime() << endl;

}

