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
#include "myfunction.h"
#include "knot_alex_table.h" // alexander polynomial table

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

//TODO 
void write_xyz(const std::string &filename, py::array_t<double> input) {
    // Get array info and validate dimensions
    py::buffer_info buf = input.request();
    
    if (buf.ndim != 3) {
        throw std::runtime_error("Input array must be 3-dimensional");
    }
    
    // Extract dimensions
    size_t nFrames = buf.shape[0];  // Number of frames
    size_t nAtoms = buf.shape[1];   // Number of atoms per frame
    size_t nCoords = buf.shape[2];  // Should be 3 (x,y,z)
    
    if (nCoords != 3) {
        throw std::runtime_error("Last dimension must be 3 (x,y,z coordinates)");
    }
    
    // Get pointer to array data
    double *ptr = static_cast<double *>(buf.ptr);
    
    // Open output file, append mode
    std::ofstream outfile(filename,std::ios::app);
    if (!outfile.is_open()) {
        throw std::runtime_error("Could not open file for writing: " + filename);
    }
    
    // Write each frame
    for (size_t frame = 0; frame < nFrames; frame++) {
        // Write number of atoms
        outfile << nAtoms << "\n";
        // Write comment line (frame number)
        outfile << "Frame " << frame + 1 << " of " << nFrames << "\n";
        
        // Write atom coordinates
        for (size_t atom = 0; atom < nAtoms; atom++) {
            // Calculate base index for current atom
            size_t base_idx = frame * (nAtoms * 3) + atom * 3;
            
            // Write atom line (assuming all atoms are carbon 'C' for simplicity)
            // Format: element_symbol x y z
            outfile << "1 "
                   << std::fixed << std::setprecision(6)
                   << ptr[base_idx] << " "      // x coordinate
                   << ptr[base_idx + 1] << " "  // y coordinate
                   << ptr[base_idx + 2] << "\n"; // z coordinate
        }
    }
    
    outfile.close();
    
    if (outfile.fail()) {
        throw std::runtime_error("Error occurred while writing file: " + filename);
    }
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

// KMT 
// 内存泄露检测
py::array_t<double> KMT_chain(py::array_t<double> input, std::string chain_type) {
    // 请求 buffer_info 以访问数组数据和维度信息
    py::buffer_info info = input.request();

    if (info.ndim != 2 || info.shape[1] != 3) {
        throw std::runtime_error("Expected a 2-dimensional array with shape (N, 3)");
    }

    ssize_t nAtoms = info.shape[0];
    ssize_t nDimension = info.shape[1];

    // 创建指向点的指针数组
    std::vector<double *> points;
    points.reserve(nAtoms);

    const double *data_ptr = static_cast<const double *>(info.ptr);

    for (ssize_t i = 0; i < nAtoms; ++i) {
        points.push_back(const_cast<double *>(data_ptr + i * nDimension));
    }

    // 根据 chain_type 调用不同的函数
    std::vector<double *> result_points;
    if (chain_type == "open") {
        result_points = KMT_open_chain(points);
    } else if (chain_type == "ring") {
        result_points = KMT(points);
    } else {
        throw std::runtime_error("Invalid chain_type. Expected 'open' or 'ring'.");
    }

    // 将结果数据从 result_points 复制到 result_data
    std::vector<double> result_data;
    for (auto &ptr : result_points) {
        for (ssize_t j = 0; j < 3; ++j) {
            result_data.push_back(ptr[j]);
        }
    }

    // 将结果转换为 NumPy 数组并返回
    std::vector<py::ssize_t> shape = {static_cast<py::ssize_t>(result_points.size()), 3};
    return py::array_t<double>(shape, result_data.data());
}

py::array_t<int> gauss_notation(py::array_t<double> input) {
    // 请求 buffer_info 以访问数组数据和维度信息
    py::buffer_info info = input.request();

    if (info.ndim != 2 || info.shape[1] != 3) {
        throw std::runtime_error("Expected a 2-dimensional array with shape (N, 3)");
    }

    ssize_t nAtoms = info.shape[0];
    ssize_t nDimension = info.shape[1];

    // 创建指向点的指针数组
    std::vector<double *> points;
    points.reserve(nAtoms);

    const double *data_ptr = static_cast<const double *>(info.ptr);

    for (ssize_t i = 0; i < nAtoms; ++i) {
        points.push_back(const_cast<double *>(data_ptr + i * nDimension));
    }

    // 调用 get_gauss_notation 函数
    std::vector<int> result_data = get_gauss_notation(points);

    // 将结果转换为 py::array_t<int> 并返回
    return py::array_t<int>(result_data.size(), result_data.data());
}

void get_alexander_map(std::string filename)
{
    std::fstream read;
    read.open(filename,std::ios::in);
    get_alexander(read);
    read.close();
}

PYBIND11_MODULE(alexander_poly, m) {
    // io part
    m.def("read_xyz", &read_xyz, "Read XYZ file and return a numpy array");
    m.def("write_xyz", &write_xyz, "Write a numpy array to a XYZ file");
    m.def("get_alexander_map", &get_alexander_map, "Get alexander map from file");
    m.def("print_alexander_map", &print_alexander_map, "Print alexander map");

    // knot type
    m.doc() = "Module for reading XYZ files with multiple frames and calculating knot type";
    m.def("calculate_knot_type", (std::vector<std::string> (*)(const std::string &)) &calculate_knot_type, "Calculate knot type from XYZ file");
    m.def("calculate_knot_type", (std::vector<std::string> (*)(py::array_t<double>)) &calculate_knot_type, "Calculate knot type from numpy array");
    
    m.def("calculate_knot_type_open_chain", (std::vector<std::string> (*)(const std::string &)) &calculate_knot_type_open_chain, "Calculate knot type from XYZ file");
    m.def("calculate_knot_type_open_chain", (std::vector<std::string> (*)(py::array_t<double>)) &calculate_knot_type_open_chain, "Calculate knot type from numpy array");
    
    // knot size
    m.def("calculate_knot_size", &calculate_knot_size, "Calculate knot size from numpy array");

    // develop
    m.def("KMT_chain", &KMT_chain, py::arg("input"), py::arg("chain_type"), "A function to compute KMT or KMT_open_chain based on chain_type");
    m.def("gauss_notation", &gauss_notation, py::arg("input"), "Calculate the Gauss notation for a given input array of points. Accepts a NumPy array of shape (N, 3) as input and returns a 1D NumPy array of integers as output.");
}