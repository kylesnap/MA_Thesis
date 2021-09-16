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

struct CellParam {
    int n;
    float b1;
};

template<typename Iterator>
void cartProduct(std::set<int> const &n, std::set<float> const &b1, Iterator out) {
    // Good, ol' fashioned nested for loops like your mom used to make (input: several vectors and an iterator
    // output: vector of tuples.)
    for (int i : n) {
        for (float j : b1) {
            CellParam temp;
            temp.n = i;
            temp.b1 = j;
            *out++ = temp;
        }
    }
}

template<typename T>
void printSet(std::set<T> set) {
    // Prints set elements
    std::for_each(set.begin(), set.end(), [](const T &e) { std::cout << e << " "; } );
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
    std::set<float> b1;
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
        n = parse_params<int>(entSizes, 1);
        b1 = parse_params<float>(entBOne, -1);
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
    char entry;
    do {
        std::cout << "Would you like to seed random? [y/n]: ";
        std::cin >> entry;
    } while (!std::cin.fail() && entry != 'y' && entry != 'n');
    (entry == 'y') ? srand(69) : srand(time(&rawtime));

    //File IO object
    do {
        std::cout << "Would you like to output the results to a CSV? [y/n]: ";
        std::cin >> entry;
    } while (!std::cin.fail() && entry != 'y' && entry != 'n');

    std::streambuf *buf;
    std::ofstream f;

    // If entry is yes, then set the stream buffer to point to an output stream to entry.
    if (entry == 'y') {
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

    // Confirm parameters
    std::cout << "~~~~~~~~~~~~~~~~~~~~ Confirm? ~~~~~~~~~~~~~~~~~~~~" << std::endl;
    std::cout << "N: ";
    printSet<int>(n);
    std::cout << "\nB1: ";
    printSet<float>(b1);

    do {
        std::cout << "\n[y/n]: ";
        std::cin >> entry;
    } while (!std::cin.fail() && entry != 'y' && entry != 'n');
    if (entry == 'n') exit(EXIT_SUCCESS);

    std::vector<CellParam> allCells;
    cartProduct(n, b1, back_inserter(allCells));

    return 0;
}