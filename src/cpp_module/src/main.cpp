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
        //std::cout<<points.size()<<std::endl;
        find_max_span(points);
        result_knottype.push_back(get_knottype_ring_faster(points));
        clean_pointer(points);
        points.clear();
    }
    read.close();

    return result_knottype;
}

std::vector<std::string> calculate_knot_type_open_chain(const std::string &filename)
{
    std::vector<std::string> result_knottype;
    std::vector<double *> points;
    std::fstream    read;
    read.open(filename,std::ios::in);

    while (read_data_cpp(points,read))
    {
        find_max_span(points);
        result_knottype.push_back(get_knottype_open_faster(points));
        clean_pointer(points);
        points.clear();
    }
    read.close();

    return result_knottype;
}
// 重构函数，输入是numpy array，输出是numpy array
std::vector<std::string> calculate_knot_type(py::array_t<double> input)
{
    // 请求buffer_info以访问数组数据和维度信息
    py::buffer_info info = input.request();

    if (info.ndim != 3) {
        throw std::runtime_error("Expected a 3-dimensional array");
    }

    // 获取指向数组数据的指针
    double *data = static_cast<double *>(info.ptr);
    ssize_t nFrames = info.shape[0];
    ssize_t nAtoms = info.shape[1];
    ssize_t nDimension = info.shape[2];

    // return value
    std::vector<std::string> result_knottype;
    std::vector<double *> points;

    for (ssize_t i = 0; i < nFrames; i++) {
        for (ssize_t j = 0; j < nAtoms; j++) {
            double *point = new double[3];
            point[0] = data[i * nAtoms * nDimension + j * nDimension + 0];
            point[1] = data[i * nAtoms * nDimension + j * nDimension + 1];
            point[2] = data[i * nAtoms * nDimension + j * nDimension + 2];
            points.push_back(point);
        }
        find_max_span(points);
        result_knottype.push_back(get_knottype_ring_faster(points));
        clean_pointer(points);
        points.clear();
    }
    return result_knottype;
}

std::vector<std::string> calculate_knot_type_open_chain(py::array_t<double> input)
{
    // 请求buffer_info以访问数组数据和维度信息
    py::buffer_info info = input.request();

    if (info.ndim != 3) {
        throw std::runtime_error("Expected a 3-dimensional array");
    }

    // 获取指向数组数据的指针
    double *data = static_cast<double *>(info.ptr);
    ssize_t nFrames = info.shape[0];
    ssize_t nAtoms = info.shape[1];
    ssize_t nDimension = info.shape[2];

    // return value
    std::vector<std::string> result_knottype;
    std::vector<double *> points;

    for (ssize_t i = 0; i < nFrames; i++) {
        for (ssize_t j = 0; j < nAtoms; j++) {
            double *point = new double[3];
            point[0] = data[i * nAtoms * nDimension + j * nDimension + 0];
            point[1] = data[i * nAtoms * nDimension + j * nDimension + 1];
            point[2] = data[i * nAtoms * nDimension + j * nDimension + 2];
            points.push_back(point);
        }
        find_max_span(points);
        result_knottype.push_back(get_knottype_open_faster(points));
        clean_pointer(points);
        points.clear();
    }
    return result_knottype;
}

// calculate knot size
std::pair<std::vector<std::string>, std::vector<std::array<int,3>>> calculate_knot_size(py::array_t<double> input,std::string chain_type)
{
    // 请求buffer_info以访问数组数据和维度信息
    py::buffer_info info = input.request();

    if (info.ndim != 3) {
        throw std::runtime_error("Expected a 3-dimensional array");
    }

    // 获取指向数组数据的指针
    double *data = static_cast<double *>(info.ptr);
    ssize_t nFrames = info.shape[0];
    ssize_t nAtoms = info.shape[1];
    ssize_t nDimension = info.shape[2];

    // return value
    std::vector<std::string> result_knottype;
    std::vector<std::array<int,3>> result_knotsize;
    std::vector<double *> points;

    for (ssize_t i = 0; i < nFrames; i++) {
        for (ssize_t j = 0; j < nAtoms; j++) {
            double *point = new double[3];
            point[0] = data[i * nAtoms * nDimension + j * nDimension + 0];
            point[1] = data[i * nAtoms * nDimension + j * nDimension + 1];
            point[2] = data[i * nAtoms * nDimension + j * nDimension + 2];
            points.push_back(point);
        }
        find_max_span(points);
        if (chain_type == "ring")
            result_knottype.push_back(get_knottype_ring_faster(points));
        else if (chain_type == "open")
            result_knottype.push_back(get_knottype_open_faster(points));

        // calculate knot size
        knot a(points);
        string temp_knottype;
        vector<int> knotSize{};

        if (chain_type == "ring")
            knotSize = a.knot_size_ring(temp_knottype,result_knottype[i]);
        else if (chain_type == "open")
            knotSize = a.knot_size(temp_knottype,result_knottype[i]);
        if (knotSize.size() <3)
            knotSize = {0,0,0};
        std::array<int,3> size = {knotSize[0], knotSize[1], knotSize[2]};
        result_knotsize.push_back(size);

        clean_pointer(points);
        points.clear();
    }
    return {result_knottype, result_knotsize};
}

PYBIND11_MODULE(alexander_poly, m) {
    m.doc() = "Module for reading XYZ files with multiple frames and calculating knot type";
    m.def("calculate_knot_type", (std::vector<std::string> (*)(const std::string &)) &calculate_knot_type, "Calculate knot type from XYZ file");
    m.def("calculate_knot_type", (std::vector<std::string> (*)(py::array_t<double>)) &calculate_knot_type, "Calculate knot type from numpy array");
    m.def("read_xyz", &read_xyz, "Read XYZ file and return a numpy array");
    m.def("calculate_knot_type_open_chain", (std::vector<std::string> (*)(const std::string &)) &calculate_knot_type_open_chain, "Calculate knot type from XYZ file");
    m.def("calculate_knot_type_open_chain", (std::vector<std::string> (*)(py::array_t<double>)) &calculate_knot_type_open_chain, "Calculate knot type from numpy array");
    m.def("calculate_knot_size", &calculate_knot_size, "Calculate knot size from numpy array");
}