// Main
// Takes in user input, prepares the list of cells to perform, then runs simulation.
// Kyle Dewsnap
// 15SEP21

#include <iostream>
#include <vector>
#include <ctime>
#include <tuple>
#include <gsl/gsl_rng.h>

#include "SimCell.h"

#define CATCH_CONFIG_RUNNER
#include "tests/catch.hpp"

#define TEST false

int run_chap4(gsl_rng *r, std::ostream *out) {
    // Runs the random factorial design of chapter four.

    int reps = 5000;
    std::tuple<float, float, float, float> trialParams = {0,1,1,0};
    std::tuple<float, float> noSat = {-1, 5};
    std::vector<float> temp;
    SimCell *cell;

    for (int i = 0; i < reps; i++) {
        temp.clear();
        int n = static_cast<int>(gsl_ran_flat(r, 250, 2500));
        std::tuple<float, float> betaM = {gsl_ran_flat(r, 0, 1), gsl_ran_flat(r, 2, 50)};
        std::tuple<float, float> betaF = {gsl_ran_flat(r, 0, 1), gsl_ran_flat(r, 2, 50)};
        std::tuple<float, float> groupProp = {gsl_ran_flat(r, 0.1, 0.9), -1};
        std::get<1>(groupProp) = 1 - std::get<0>(groupProp);

        // 25% change of no satisficing in each group (~% of no satisficing in both groups)
        betaM = gsl_rng_uniform(r) < 0.25 ? noSat : betaM;
        betaF = gsl_rng_uniform(r) < 0.35 ? noSat : betaF;

        cell = new SimCell(n, 1, betaM, betaF, groupProp, trialParams, r, false);
        cell->toVec(temp, true);
        *out << cell->run();
    }

    return 0;
}

int run_chap5(gsl_rng *r, std::ostream *out) {
    // Runs the blocked random factorial design of chapter five.

    int sampleSize[3] = {250, 1000, 2500};
    float betaModes[3] = {0, 0.1, 0.25};
    float betaK[3] = {10, 25, 50};

    // Make a list of beta distribution params
    std::vector<std::tuple<float, float>> betaDists;
    for (float i : betaModes) {
        for (float j : betaK) {
            betaDists.emplace_back(std::make_tuple(i, j));
        }
    }

    int reps = 100;
    std::tuple<float, float, float, float> trialParams = {0,1,0,1};
    std::vector<float> temp;
    SimCell *cell;

    for (int n : sampleSize) {
        for (auto betaM : betaDists) {
            for (auto betaF : betaDists) {
                for (int i = 0; i < reps; i++) {
                    // Get two numbers, each bigger than 0.1, smaller than 0.9.
                    float propMale = 1, propFemale = 1;
                    while (propMale + propFemale >= 1) {
                        propMale = (float) gsl_ran_flat(r, 0.1, 0.9);
                        propFemale = (float) gsl_ran_flat(r, 0.1, 0.9);
                    }

                    // Run simulation.
                    temp.clear();
                    std::tuple<float, float> groupProp = {propMale, propFemale};
                    cell = new SimCell(n, 1, betaM, betaF, groupProp, trialParams, r, false);
                    cell->toVec(temp, true);
                    *out << cell->run();
                }
            }
        }

    }

    return 0;
}

int main() {
    char entry;

    std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
    std::cout << "Welcome to KD's Simulation of Linear Regression!" << std::endl;

    // This macro will run tests if called from Catch.hpp (Macro set earlier)
    if (TEST) {
        return Catch::Session().run();
    }

    // Prepare random and ask for seeding.
    gsl_rng_env_setup();
    const gsl_rng_type *T = gsl_rng_default;
    gsl_rng *r = gsl_rng_alloc(T);

    do {
        std::cout << "Seed random? [y/n]: ";
        std::cin >> entry;
    } while (!std::cin.fail() && entry != 'y' && entry != 'n');
    if (entry == 'y')  gsl_rng_set(r, 69);

    // File IO Prep.
    do {
        std::cout << "Print to File? [y/n]: ";
        std::cin >> entry;
    } while (!std::cin.fail() && entry != 'y' && entry != 'n');

    std::ostream *out;
    if (entry == 'y') {
        time_t rawTime;
        struct tm *timeInfo;
        time(&rawTime);
        timeInfo = localtime(&rawTime);

        char fChar[32];
        strftime(fChar, 32, "simOut_%d%h%y_%H%M%S.csv", timeInfo);
        std::cout << "File will be written in cur. dir: " << fChar << std::endl;

        out = new std::ofstream(fChar);
    } else {
        out = &std::cout;
    }

    // Run either sim for chap 4 or chap 5
    do {
        std::cout << "Which chapter would you like to run? [4/5]: ";
        std::cin >> entry;
    } while (!std::cin.fail() && entry != '4' && entry != '5');

    *out << "I,N,NM,NF,NX,BM_A,BM_B,BF_A,BF_B,"
            "BTRUE_0,BTRUE_1,BTRUE_F,BTRUE_X,ERR_VAR,"
            "BHAT_0,BHAT_1,BHAT_F,BSE_0,BSE_1,BSE_F,RSQ,"
            "NBHAT_0,NBHAT_1,NBHAT_F,NBSE_0,NBSE_1,NBSE_F,NRSQ,NS,OBS_F" << std::endl;

    if (entry == '4') run_chap4(r, out);
    if (entry == '5') run_chap5(r, out);
    // CHAP 5 LINE HERE.

    if (out != &std::cout) delete out;
    return 0;
}
