/*
    Pybind11 function, input: filename, output: numpy array of knot type
*/

#include <vector>
#include <string>
#include <iostream>
#include <algorithm>
#include <fstream>
#include <map>
#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <unistd.h>
#include <getopt.h>
#include "knot.h"
#include "knottype.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
namespace py = pybind11;

// global variable
int flag_ring_open=1; //1 for ring

struct Atom {
    std::string type;
    double x, y, z;
};

struct Frame {
    std::vector<Atom> atoms;
};

py::array_t<double> atoms_to_numpy(const std::vector<Frame>& frames) {
    size_t max_atoms = 0;
    for (const auto& frame : frames) {
        max_atoms = std::max(max_atoms, frame.atoms.size());
    }

    size_t num_frames = frames.size();
    std::vector<double> data(num_frames * max_atoms * 3, std::numeric_limits<double>::quiet_NaN()); // Fill with NaN

    for (size_t i = 0; i < num_frames; i++) {
        const auto& atoms = frames[i].atoms;
        for (size_t j = 0; j < atoms.size(); j++) {
            size_t index = i * max_atoms * 3 + j * 3;
            data[index + 0] = atoms[j].x;
            data[index + 1] = atoms[j].y;
            data[index + 2] = atoms[j].z;
        }
    }
    std::vector<py::ssize_t> shape = {num_frames, max_atoms, 3};
    return py::array_t<double>(shape, data.data());
}

py::array_t<double> read_xyz(const std::string &filename) {
    std::vector<Frame> frames;
    std::ifstream file(filename);
    std::string line;
    int num_atoms;

    while (file >> num_atoms) {
        std::getline(file, line);  // Finish reading the num_atoms line
        std::getline(file, line);  // Skip the comment line
        Frame frame;
        for (int i = 0; i < num_atoms; i++) {
            std::getline(file, line);
            std::istringstream iss(line);
            Atom atom;
            if (!(iss >> atom.type >> atom.x >> atom.y >> atom.z)) {
                break; // Error in parsing line
            }
            frame.atoms.push_back(atom);
        }
        frames.push_back(frame);
    }
    return atoms_to_numpy(frames);
}

std::vector<std::string> calculate_knot_type(const std::string &filename) 
{
    std::vector<std::string> result_knottype;
    std::vector<double *> points;
    std::fstream    read;
    read.open(filename,std::ios::in);

    //print_alexander_map();
    
    while(read_data_cpp(points,read))
    {
        std::cout<<points.size()<<std::endl;
        find_max_span(points);
        result_knottype.push_back(get_knottype_ring_faster(points));
        clean_pointer(points);
        points.clear();
    }
    read.close();

    return result_knottype;
}

PYBIND11_MODULE(alexander_poly, m) {
    m.doc() = "Module for reading XYZ files with multiple frames and calculating knot type";
    m.def("calculate_knot_type", &calculate_knot_type, "Calculate knot type from a file");
    m.def("read_xyz", &read_xyz, "Read XYZ file and return a numpy array");
}