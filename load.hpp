// load.hpp
#ifndef LOAD_HPP
#define LOAD_HPP

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>

std::vector<double> read_csv_1d(const std::string& filename) {
    std::vector<double> data;
    std::ifstream file(filename);
    std::string line;

    if (!file.is_open()) {
        std::cerr << "Failed to open file: " << filename << "\n";
        return data;
    }

    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string value;
        while (std::getline(ss, value, ',')) {
            data.push_back(std::stod(value));
        }
    }

    return data;
}

#endif 