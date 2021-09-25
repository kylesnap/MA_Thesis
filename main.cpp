// Main
// Takes in user input, prepares the list of cells to perform, then runs simulation.
// Kyle Dewsnap
// 15SEP21

#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <set>
#include <tuple>

#include "SimCell.h"

template<typename Iterator>
void cartProduct(std::set<int> const &n, std::set<float> const &b1, Iterator out) {
    // Good, ol' fashioned nested for loops like your mom used to make (input: several vectors and an iterator
    // output: vector of tuples.)
    for (int i: n) {
        for (float j: b1) {
            CellParam temp;
            temp.n = i;
            temp.testFloat = j;
            *out++ = temp;
        }
    }
}

template<typename T>
void printSet(std::set<T> set) {
    // Prints set elements
    std::for_each(set.begin(), set.end(), [](const T &e) { std::cout << e << " "; });
}

template<typename T>
T strToType(std::string const &str, int range) {
    // takes in string, transforms to type 'T' and checks whether range is met.

    T value = (typeid(T) == typeid(int)) ? std::stoi(str) :
              (typeid(T) == typeid(float)) ? std::stof(str) :
              throw std::bad_typeid();

    if (value > 0 && range > 0 || value >= 0 && range == 0 || range < 0) return value;
    throw std::out_of_range("One or more arguments did not follow the range specifications.");
}

template<typename T>
std::set<T> parse_params(std::string str, int range) { // Here
    // Takes in user entered parameters, cleans them, then outputs 2D array of numeric parameters.
    std::vector<T> entry;

    while (!empty(str)) {
        std::size_t found = str.find(',');
        entry.push_back(strToType<T>(str.substr(0, found), range));
        str = found == std::string::npos ? "" : str.substr(++found);
    }

    return std::set<T>(entry.begin(), entry.end());
}

int main() {
    std::set<int> n;
    std::set<float> testFloat;
    std::string entSizes, entFloat;

    std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
    std::cout << "Welcome to KD's Simulation of Linear Regression!" << std::endl;

    // Ask whether parameters should be acquired, or if this is a test run.
    char entry;
    do {
        std::cout << "Is this a test? (i.e., use dummy parameters) [y/n]: ";
        std::cin >> entry;
    } while (!std::cin.fail() && entry != 'y' && entry != 'n');

    if (entry == 'y') {
        std::cout << "Skipping input steps..." << std::endl;
        n = {5, 10};
        testFloat = {1, 2};
    } else {
        // Asking and processing simulation parameters
        std::cout << "Please enter the following parameters separated by commas:" << std::endl;
        std::cout << "Sample Sizes (Must be larger than zero): ";
        getline(std::cin, entSizes);
        std::cout << "Test Parameter of Floats: ";
        std::getline(std::cin, entFloat);
        try {
            n = parse_params<int>(entSizes, 1);
            testFloat = parse_params<float>(entFloat, -1);
        }
        catch (std::invalid_argument &e) {
            std::cout << "Error in parsing your arguments: "
                         "Please only use numerals and commas in input line." << std::endl;
            main();
        }
        catch (std::out_of_range &e) {
            std::cout << "Error in parsing your arguments: Please follow the range guidelines." << std::endl;
            main();
        }
        catch (...) {
            std::cout << "Unspecified error occurred. Check typing and weep." << std::endl;
            exit(EXIT_FAILURE);
        }
    }

    time_t rawTime;
    struct tm *timeInfo;

    time(&rawTime);
    timeInfo = localtime(&rawTime);

    // Ask whether random should be seeded (for testing!)
    // TODO: Rework this function to use the GSL random generator
    do {
        std::cout << "Would you like to seed random? [y/n]: ";
        std::cin >> entry;
    } while (!std::cin.fail() && entry != 'y' && entry != 'n');
    (entry == 'y') ? srand(69) : srand(time(&rawTime));

    //File IO object
    do {
        std::cout << "Would you like to output the results to a CSV? [y/n]: ";
        std::cin >> entry;
    } while (!std::cin.fail() && entry != 'y' && entry != 'n');

    std::string fName;
    // If entry is yes, then set the stream buffer to point to an output stream to entry.
    if (entry == 'y') {
        char fChar[32];
        strftime(fChar, 32, "simOut_%d%h%y_%H%M%S.csv", timeInfo);
        fName = fChar;
        std::cout << "File will be written in cur. dir: " << fName << std::endl;
    } else {
        // Else, set the stream to just point to cout!
        fName[0] = 0;
    }

    // Confirm parameters
    std::cout << "~~~~~~~~~~~~~~~~~~~~ Confirm? ~~~~~~~~~~~~~~~~~~~~" << std::endl;
    std::cout << "N: ";
    printSet<int>(n);
    std::cout << "\nTest Float: ";
    printSet<float>(testFloat);

    do {
        std::cout << "\n[y/n]: ";
        std::cin >> entry;
    } while (!std::cin.fail() && entry != 'y' && entry != 'n');
    if (entry == 'n') exit(EXIT_SUCCESS);

    std::vector<CellParam> allCells;
    cartProduct(n, testFloat, back_inserter(allCells));

    SimCell *currCell;
    for (auto i : allCells) {
        try {
            currCell = new SimCell(i, fName);
            currCell->run();
        }
        catch (std::bad_function_call &e) {
            std::cout << "Cell failed to construct. Skipping.";
            continue;
        }
    }

    return 0;
}