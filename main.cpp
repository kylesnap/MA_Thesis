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

    std::vector<int> trialN = {100, 500, 1000};
    std::vector<float> trialVar = {1, 2, 3};
    std::vector<std::tuple<float, float>> trialBetaP = {{-1,50}};
    std::vector<std::tuple<float, float>> trialBetaQ = {{-1,50}};
    std::vector<std::tuple<float, float>> trialPropG = {{0.5,0.5}};
    std::vector<std::tuple<float, float, float, float>> trialParams = {{1,0,0,0}};

    *out << "I,N,NP,NQ,NX,BP_A,BP_B,BQ_A,BQ_B,"
           "BTRUE_0,BTRUE_1,BTRUE_Q,BTRUE_X,ERR_VAR,"
           "BHAT_0,BHAT_1,BHAT_Q,BSE_0,BSE_1,BSE_Q,RSQ" << std::endl;

    // God forgive me for this loop
    SimCell *cell;
    std::vector<float> ran;
    for (auto n : trialN) {
        for (auto var : trialVar) {
            for (auto bP : trialBetaP) {
                for (auto bQ : trialBetaQ) {
                    for (auto pG : trialPropG) {
                        for (auto prm : trialParams) {
                            cell = new SimCell(n, var, bP, bQ, pG, prm, r);
                            cell->toVec(ran, true);
                            *out << cell->run() << std::endl;
                        }
                    }
                }
            }
        }
    }

    if (out != &std::cout) delete out;
    return 0;
}