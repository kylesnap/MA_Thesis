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

template<typename Iterator>
void cartProduct(std::set<int> const &n, std::set<float> const &b1, Iterator out) {
    // Good, ol' fashioned nested for loops like your mom used to make (input: several vectors and an iterator
    // output: vector of tuples.)
    for (int i : n) {
        for (float j : b1) {
            *out++ = std::make_tuple(i, j);
        }
    }
}

template<typename T>
void strToType(std::string const &str, T *value, int range) {
    // takes in string, transforms to type 'T' and checks whether range is met.

    if (typeid(T) == typeid(int)) {
        *value = std::stoi(str);
    } else if (typeid(T) == typeid(float)) {
        *value = std::stof(str);
    } else {
        throw std::bad_typeid();
    }

    if (*value > 0 && range > 0 || *value >= 0 && range == 0 || range < 0) return;
    throw std::out_of_range("One or more arguments did not follow the range specifications.");
}

template<typename T>
void parse_params(std::string str, std::vector<T> *result,  int range) { // Here
    // Takes in user entered parameters, cleans them, then outputs 2D array of numeric parameters.
    T temp;

    while (!empty(str)) {
        std::size_t found = str.find(',');
        if (found != std::string::npos) {
            strToType<T>(str.substr(0, found), &temp, range);
            result->push_back(temp);
            str = str.substr(++found);
        } else {
            strToType<T>(str.substr(0, found), &temp, range);
            result->push_back(temp);
            str = "";
        }
    }
}

int main() {
    std::vector<int> sampleSizes;
    std::vector<float> betaOne;
    std::string entSizes, entBOne;

    std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
    std::cout << "Welcome to KD's Simulation of Logistic Regression!" << std::endl;

    // Asking and processing simulation parameters
    std::cout << "Please enter the following parameters separated by commas:" << std::endl;
    std::cout << "Sample Sizes (Must be larger than zero): ";
    getline (std::cin, entSizes);
    std::cout << "Beta 1 (First Parameter): ";
    getline (std::cin, entBOne);
    try {
        parse_params<int>(entSizes, &sampleSizes, 1); // The address of this thing is sent to...
        parse_params<float>(entBOne, &betaOne, -1);
    }
    catch (std::invalid_argument& e) {
        std::cout << "Error in parsing your arguments: Please only use numerals and commas in input line." << std::endl;
        main();
    }
    catch (std::out_of_range& e) {
        std::cout << "Error in parsing your arguments: Please follow the range guidelines." << std::endl;
        main();
    }
    catch (...) {
        std::cout << "Unspecified error occurred. Check typing and weep." << std::endl;
        exit(EXIT_FAILURE);
    }

    time_t rawtime;
    struct tm * timeinfo;

    time (&rawtime);
    timeinfo = localtime(&rawtime);

    // Ask whether random should be seeded (for testing!)
    char seed;
    do {
        std::cout << "Would you like to seed random? [y/n]: ";
        std::cin >> seed;
    } while (!std::cin.fail() && seed != 'y' && seed != 'n');
    (seed == 'y') ? srand(69) : srand(time(&rawtime));

    //File IO object
    char file;
    do {
        std::cout << "Would you like to output the results to a CSV? [y/n]: ";
        std::cin >> file;
    } while (!std::cin.fail() && file != 'y' && file != 'n');

    std::streambuf *buf;
    std::ofstream f;

    // If file is yes, then set the stream buffer to point to an output stream to file.
    if (file == 'y') {
        char fName [32];
        strftime(fName, 32, "simOut_%d%h%y_%H%M%S.csv", timeinfo);
        f.open(fName);
        buf = f.rdbuf();
        std::cout << "File will be written in cur. dir: " << fName << std::endl;
    } else {
        // Else, set the stream to just point to cout!
        buf = std::cout.rdbuf();
    }

    // Finally, refer the stream buffer to this output object.
    std::ostream out(buf);
    out << "Test!" << std::endl;

    // Makes table of inputs!
    std::set<int> n = std::set<int>(sampleSizes.begin(), sampleSizes.end());
    std::set<float> b1 = std::set<float>(betaOne.begin(), betaOne.end());

    std::vector<std::tuple<int, float>> cellParams;
    cartProduct(n, b1, back_inserter(cellParams));

    return 0;
}