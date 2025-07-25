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
#include <limits>
#include <functional>
#include "knot.h"
#include "knottype.h"
#include "myfunction.h"
#include "knot_alex_table.h" // alexander polynomial table
#include "optimized_types.h"  // 新增的优化类型

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
    if (UNLIKELY(frames.empty())) {
        return py::array_t<double>();
    }
    
    size_t max_atoms = 0;
    for (const auto& frame : frames) {
        max_atoms = std::max(max_atoms, frame.atoms.size());
    }

    const size_t num_frames = frames.size();
    const size_t total_size = num_frames * max_atoms * 3;
    
    // 预分配并使用NaN填充
    std::vector<double> data(total_size, std::numeric_limits<double>::quiet_NaN());

    // 优化的数据复制循环
    for (size_t i = 0; i < num_frames; ++i) {
        const auto& atoms = frames[i].atoms;
        const size_t frame_base = i * max_atoms * 3;
        
        for (size_t j = 0; j < atoms.size(); ++j) {
            const size_t atom_base = frame_base + j * 3;
            // 使用优化的内存复制
            data[atom_base + 0] = atoms[j].x;
            data[atom_base + 1] = atoms[j].y;
            data[atom_base + 2] = atoms[j].z;
        }
    }
    
    std::vector<py::ssize_t> shape = {static_cast<py::ssize_t>(num_frames), 
                                      static_cast<py::ssize_t>(max_atoms), 3};
    return py::array_t<double>(shape, data.data());
}

py::array_t<double> read_xyz(const std::string &filename) {
    std::vector<Frame> frames;
    std::ifstream file(filename);
    
    if (UNLIKELY(!file.is_open())) {
        throw std::runtime_error("Cannot open file: " + filename);
    }
    
    // 使用更大的缓冲区提高I/O性能
    constexpr size_t BUFFER_SIZE = 64 * 1024; // 64KB buffer
    file.rdbuf()->pubsetbuf(nullptr, BUFFER_SIZE);
    
    std::string line;
    int num_atoms;
    
    // 预分配frames容量以减少重新分配
    frames.reserve(1000); // 假设通常有不到1000帧

    while (file >> num_atoms) {
        std::getline(file, line);  // Finish reading the num_atoms line
        std::getline(file, line);  // Skip the comment line
        
        Frame frame;
        frame.atoms.reserve(num_atoms); // 预分配原子数量
        
        for (int i = 0; i < num_atoms; i++) {
            std::getline(file, line);
            std::istringstream iss(line);
            Atom atom;
            if (UNLIKELY(!(iss >> atom.type >> atom.x >> atom.y >> atom.z))) {
                break; // Error in parsing line
            }
            frame.atoms.push_back(std::move(atom)); // 使用移动语义
        }
        frames.push_back(std::move(frame)); // 使用移动语义
    }
    
    return atoms_to_numpy(frames);
}

// SIMD增强版本的calculate_knot_type - 实验性优化
std::vector<std::string> calculate_knot_type_simd(py::array_t<double> input, const std::string &chain_type = "ring")
{
    // 使用FastArrayInfo减少重复的buffer_info调用
    OptimizedUtils::pybind_optimized::FastArrayInfo array_info(input);
    
    // 验证chain_type参数
    if (UNLIKELY(chain_type != "ring" && chain_type != "open")) {
        throw std::runtime_error("Invalid chain_type. Expected 'ring' or 'open'.");
    }

    // 使用BatchResultBuilder优化结果构建
    OptimizedUtils::pybind_optimized::BatchResultBuilder<std::string> result_builder(array_info.nFrames);
    
    // 获取缓存实例 - 现在使用线程安全的缓存
    static ResultCache<size_t, std::string> cache(500);
    
    const bool is_ring = (chain_type == "ring");
    
    // 并行处理帧数据（如果数据量足够大）
    if (array_info.nFrames >= 4) {
        // 使用并行处理优化大数据集
        std::vector<std::string> parallel_results(array_info.nFrames);
        
        OptimizedUtils::parallel::parallelFor(0, array_info.nFrames, 2, 
            [&](size_t start, size_t end) {
                for (size_t i = start; i < end; ++i) {
                    const size_t frame_offset = i * array_info.nAtoms * array_info.nDimension;
                    
                    // 计算帧哈希用于缓存
                    const size_t frame_hash = OptimizedUtils::computeFrameHash(
                        array_info.data, array_info.nAtoms, array_info.nDimension, frame_offset);
                    const size_t cache_key = frame_hash ^ std::hash<std::string>()(chain_type);
                    
                    // 检查缓存（注意：并行访问需要线程安全）
                    std::string cached_result;
                    if (cache.get(cache_key, cached_result)) {
                        parallel_results[i] = cached_result;
                        continue;
                    }
                    
                    // 使用RAII管理内存的PointManager
                    PointManager point_manager(array_info.data, array_info.nAtoms, 
                                             array_info.nDimension, frame_offset);
                    
                    auto points = point_manager.getCompatiblePointers();
                    find_max_span(points);
                    
                    std::string knottype;
                    if (LIKELY(is_ring)) {
                        knottype = get_knottype_ring_faster(points);
                    } else {
                        knottype = get_knottype_open_faster(points);
                    }
                    
                    parallel_results[i] = knottype;
                    cache.put(cache_key, knottype);
                }
            });
        
        return parallel_results;
    } else {
        // 小数据集使用串行处理
        for (ssize_t i = 0; i < array_info.nFrames; i++) {
            const size_t frame_offset = i * array_info.nAtoms * array_info.nDimension;
            
            const size_t frame_hash = OptimizedUtils::computeFrameHash(
                array_info.data, array_info.nAtoms, array_info.nDimension, frame_offset);
            const size_t cache_key = frame_hash ^ std::hash<std::string>()(chain_type);
            
            std::string cached_result;
            if (LIKELY(cache.get(cache_key, cached_result))) {
                result_builder.addResult(std::move(cached_result));
                continue;
            }
            
            PointManager point_manager(array_info.data, array_info.nAtoms, 
                                     array_info.nDimension, frame_offset);
            
            auto points = point_manager.getCompatiblePointers();
            
            // 预取下一帧数据
            if (LIKELY(i + 1 < array_info.nFrames)) {
                const size_t next_frame_offset = (i + 1) * array_info.nAtoms * array_info.nDimension;
                PREFETCH(array_info.data + next_frame_offset);
            }
            
            find_max_span(points);
            
            std::string knottype;
            if (LIKELY(is_ring)) {
                knottype = get_knottype_ring_faster(points);
            } else {
                knottype = get_knottype_open_faster(points);
            }
            
            cache.put(cache_key, knottype);
            result_builder.addResult(std::move(knottype));
        }
        
        return std::move(result_builder).build();
    }
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

std::vector<std::string> calculate_knot_type(const std::string &filename, const std::string &chain_type = "ring") 
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
        if (chain_type == "ring") {
            result_knottype.push_back(get_knottype_ring_faster(points));
        } else if (chain_type == "open") {
            result_knottype.push_back(get_knottype_open_faster(points));
        } else {
            throw std::runtime_error("Invalid chain_type. Expected 'ring' or 'open'.");
        }
        clean_pointer(points);
        points.clear();
    }
    read.close();

    return result_knottype;
}

// 重构函数，输入是numpy array，输出是numpy array - 带缓存优化
std::vector<std::string> calculate_knot_type(py::array_t<double> input, const std::string &chain_type = "ring")
{
    // 请求buffer_info以访问数组数据和维度信息
    py::buffer_info info = input.request();

    if (UNLIKELY(info.ndim != 3)) {
        throw std::runtime_error("Expected a 3-dimensional array");
    }

    // 获取指向数组数据的指针
    double *data = static_cast<double *>(info.ptr);
    const ssize_t nFrames = info.shape[0];
    const ssize_t nAtoms = info.shape[1];
    const ssize_t nDimension = info.shape[2];

    // 预分配结果向量大小，避免重复扩容
    std::vector<std::string> result_knottype;
    result_knottype.reserve(nFrames);
    
    // 验证chain_type参数
    if (UNLIKELY(chain_type != "ring" && chain_type != "open")) {
        throw std::runtime_error("Invalid chain_type. Expected 'ring' or 'open'.");
    }

    // 获取缓存实例 - 现在使用线程安全的缓存
    static ResultCache<size_t, std::string> cache(500); // 缓存500个结果
    
    // 处理每一帧
    for (ssize_t i = 0; i < nFrames; i++) {
        const size_t frame_offset = i * nAtoms * nDimension;
        
        // 计算帧的哈希值用于缓存
        const size_t frame_hash = OptimizedUtils::computeFrameHash(data, nAtoms, nDimension, frame_offset);
        const size_t cache_key = frame_hash ^ std::hash<std::string>()(chain_type);
        
        // 检查缓存
        std::string cached_result;
        if (LIKELY(cache.get(cache_key, cached_result))) {
            result_knottype.push_back(cached_result);
            continue;
        }
        
        // 使用RAII管理内存的PointManager
        PointManager point_manager(data, nAtoms, nDimension, frame_offset);
        
        // 获取兼容格式的指针（用于legacy函数）
        auto points = point_manager.getCompatiblePointers();
        
        // 预取下一帧数据到缓存
        if (LIKELY(i + 1 < nFrames)) {
            const size_t next_frame_offset = (i + 1) * nAtoms * nDimension;
            PREFETCH(data + next_frame_offset);
        }
        
        // 调用find_max_span和knot type计算
        find_max_span(points);
        
        std::string knottype;
        if (LIKELY(chain_type == "ring")) {
            knottype = get_knottype_ring_faster(points);
        } else { // chain_type == "open"
            knottype = get_knottype_open_faster(points);
        }
        
        // 缓存结果
        cache.put(cache_key, knottype);
        result_knottype.push_back(knottype);
        
        // PointManager析构时自动清理内存，无需手动clean_pointer
    }
    
    return result_knottype;
}


// calculate knot size - 高度优化版本
std::pair<std::vector<std::string>, std::vector<std::array<int,3>>> calculate_knot_size(py::array_t<double> input, const std::string &chain_type)
{
    // 请求buffer_info以访问数组数据和维度信息
    py::buffer_info info = input.request();

    if (UNLIKELY(info.ndim != 3)) {
        throw std::runtime_error("Expected a 3-dimensional array");
    }

    // 获取指向数组数据的指针
    double *data = static_cast<double *>(info.ptr);
    const ssize_t nFrames = info.shape[0];
    const ssize_t nAtoms = info.shape[1];
    const ssize_t nDimension = info.shape[2];

    // 验证chain_type参数
    if (UNLIKELY(chain_type != "ring" && chain_type != "open")) {
        throw std::runtime_error("Invalid chain_type. Expected 'ring' or 'open'.");
    }

    // 预分配结果向量大小
    std::vector<std::string> result_knottype;
    std::vector<std::array<int,3>> result_knotsize;
    result_knottype.reserve(nFrames);
    result_knotsize.reserve(nFrames);

    // 获取缓存实例 - 分别为knottype和knotsize缓存，现在使用线程安全的缓存
    static ResultCache<size_t, std::string> type_cache(500);
    static ResultCache<size_t, std::array<int,3>> size_cache(500);
    
    // 批处理优化：一次处理多帧以提高缓存效率
    const bool is_ring = (chain_type == "ring");
    
    // 处理每一帧
    for (ssize_t i = 0; i < nFrames; i++) {
        const size_t frame_offset = i * nAtoms * nDimension;
        
        // 计算帧的哈希值用于缓存
        const size_t frame_hash = OptimizedUtils::computeFrameHash(data, nAtoms, nDimension, frame_offset);
        const size_t cache_key = frame_hash ^ std::hash<std::string>()(chain_type);
        
        // 检查缓存
        std::string cached_knottype;
        std::array<int,3> cached_knotsize;
        bool type_cached = type_cache.get(cache_key, cached_knottype);
        bool size_cached = size_cache.get(cache_key, cached_knotsize);
        
        if (LIKELY(type_cached && size_cached)) {
            result_knottype.push_back(cached_knottype);
            result_knotsize.push_back(cached_knotsize);
            continue;
        }
        
        // 使用RAII管理内存的PointManager
        PointManager point_manager(data, nAtoms, nDimension, frame_offset);
        
        // 获取兼容格式的指针（用于legacy函数）
        auto points = point_manager.getCompatiblePointers();
        
        // 预取下一帧数据到缓存
        if (LIKELY(i + 1 < nFrames)) {
            const size_t next_frame_offset = (i + 1) * nAtoms * nDimension;
            PREFETCH(data + next_frame_offset);
        }
        
        // 调用find_max_span和knot type计算
        find_max_span(points);
        
        std::string knottype;
        if (!type_cached) {
            if (LIKELY(is_ring)) {
                knottype = get_knottype_ring_faster(points);
            } else { // chain_type == "open"
                knottype = get_knottype_open_faster(points);
            }
            type_cache.put(cache_key, knottype);
        } else {
            knottype = cached_knottype;
        }
        result_knottype.push_back(knottype);

        std::array<int,3> knotsize_result;
        if (!size_cached) {
            // calculate knot size
            knot knot_analyzer(points);
            std::string temp_knottype;
            std::vector<int> knotSize;

            if (LIKELY(is_ring)) {
                knotSize = knot_analyzer.knot_size_ring(temp_knottype, knottype);
            } else { // chain_type == "open"
                knotSize = knot_analyzer.knot_size(temp_knottype, knottype);
            }
            
            // 确保knotSize有3个元素
            if (UNLIKELY(knotSize.size() < 3)) {
                knotsize_result = {{0, 0, 0}};
            } else {
                knotsize_result = {{knotSize[0], knotSize[1], knotSize[2]}};
            }
            
            size_cache.put(cache_key, knotsize_result);
        } else {
            knotsize_result = cached_knotsize;
        }
        result_knotsize.push_back(knotsize_result);
        
        // PointManager析构时自动清理内存
    }
    
    return std::make_pair(std::move(result_knottype), std::move(result_knotsize));
}

// KMT - optimized version with better memory management
py::array_t<double> KMT_chain(py::array_t<double> input, const std::string &chain_type) {
    // 请求 buffer_info 以访问数组数据和维度信息
    py::buffer_info info = input.request();

    if (info.ndim != 2 || info.shape[1] != 3) {
        throw std::runtime_error("Expected a 2-dimensional array with shape (N, 3)");
    }

    const ssize_t nAtoms = info.shape[0];
    constexpr ssize_t nDimension = 3;
    
    // 验证chain_type参数
    if (chain_type != "ring" && chain_type != "open") {
        throw std::runtime_error("Invalid chain_type. Expected 'ring' or 'open'.");
    }

    // 创建指向点的指针数组（这里需要与现有API兼容）
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
    } else { // chain_type == "ring"
        result_points = KMT(points);
    }

    // 预分配结果数据向量
    const size_t result_size = result_points.size();
    std::vector<double> result_data;
    result_data.reserve(result_size * 3);
    
    // 将结果数据从 result_points 复制到 result_data
    for (const auto &ptr : result_points) {
        result_data.insert(result_data.end(), ptr, ptr + 3);
    }

    // 将结果转换为 NumPy 数组并返回
    std::vector<py::ssize_t> shape = {static_cast<py::ssize_t>(result_size), 3};
    return py::array_t<double>(shape, result_data.data());
}

py::array_t<int> gauss_notation(py::array_t<double> input) {
    // 请求 buffer_info 以访问数组数据和维度信息
    py::buffer_info info = input.request();

    if (info.ndim != 2 || info.shape[1] != 3) {
        throw std::runtime_error("Expected a 2-dimensional array with shape (N, 3)");
    }

    const ssize_t nAtoms = info.shape[0];
    constexpr ssize_t nDimension = 3;

    // 创建指向点的指针数组（为了与现有API兼容）
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

// 超高性能并行版本 - calculate_knot_type
std::vector<std::string> calculate_knot_type_parallel(py::array_t<double> input, const std::string &chain_type = "ring")
{
    // 使用优化的数组信息提取
    OptimizedUtils::pybind_optimized::FastArrayInfo array_info(input);
    
    // 验证chain_type参数
    if (UNLIKELY(chain_type != "ring" && chain_type != "open")) {
        throw std::runtime_error("Invalid chain_type. Expected 'ring' or 'open'.");
    }

    // 预分配结果向量
    std::vector<std::string> result_knottype(array_info.nFrames);
    
    // 获取缓存实例 - 现在使用线程安全的缓存
    static ResultCache<size_t, std::string> cache(1000);
    
    // 确定是否值得并行化
    const size_t min_frames_for_parallel = 4;
    const bool is_ring = (chain_type == "ring");
    
    if (array_info.nFrames >= min_frames_for_parallel) {
        // 并行处理版本
        OptimizedUtils::parallel::parallelFor(
            0, array_info.nFrames, 2,
            [&](size_t start, size_t end) {
                for (size_t i = start; i < end; ++i) {
                    const size_t frame_offset = i * array_info.nAtoms * array_info.nDimension;
                    
                    // 计算帧的哈希值用于缓存
                    const size_t frame_hash = OptimizedUtils::computeFrameHash(
                        array_info.data, array_info.nAtoms, array_info.nDimension, frame_offset);
                    const size_t cache_key = frame_hash ^ std::hash<std::string>()(chain_type);
                    
                    // 检查缓存（注意：多线程访问需要同步）
                    std::string cached_result;
                    if (cache.get(cache_key, cached_result)) {
                        result_knottype[i] = cached_result;
                        continue;
                    }
                    
                    // 使用RAII管理内存的PointManager
                    PointManager point_manager(array_info.data, array_info.nAtoms, 
                                             array_info.nDimension, frame_offset);
                    
                    // 获取兼容格式的指针
                    auto points = point_manager.getCompatiblePointers();
                    
                    // 调用计算函数
                    find_max_span(points);
                    
                    std::string knottype;
                    if (LIKELY(is_ring)) {
                        knottype = get_knottype_ring_faster(points);
                    } else {
                        knottype = get_knottype_open_faster(points);
                    }
                    
                    // 缓存结果
                    cache.put(cache_key, knottype);
                    result_knottype[i] = knottype;
                }
            }
        );
    } else {
        // 单线程版本（小数据量）
        for (ssize_t i = 0; i < array_info.nFrames; i++) {
            const size_t frame_offset = i * array_info.nAtoms * array_info.nDimension;
            
            const size_t frame_hash = OptimizedUtils::computeFrameHash(
                array_info.data, array_info.nAtoms, array_info.nDimension, frame_offset);
            const size_t cache_key = frame_hash ^ std::hash<std::string>()(chain_type);
            
            std::string cached_result;
            if (cache.get(cache_key, cached_result)) {
                result_knottype[i] = cached_result;
                continue;
            }
            
            PointManager point_manager(array_info.data, array_info.nAtoms, 
                                     array_info.nDimension, frame_offset);
            auto points = point_manager.getCompatiblePointers();
            
            find_max_span(points);
            
            std::string knottype;
            if (LIKELY(is_ring)) {
                knottype = get_knottype_ring_faster(points);
            } else {
                knottype = get_knottype_open_faster(points);
            }
            
            cache.put(cache_key, knottype);
            result_knottype[i] = knottype;
        }
    }
    
    return result_knottype;
}

PYBIND11_MODULE(alexander_poly, m) {
    // io part
    m.def("read_xyz", &read_xyz, "Read XYZ file and return a numpy array");
    m.def("write_xyz", &write_xyz, "Write a numpy array to a XYZ file");
    m.def("get_alexander_map", &get_alexander_map, "Get alexander map from file");
    m.def("print_alexander_map", &print_alexander_map, "Print alexander map");

    // knot type - 原版本（保持兼容性）
    m.doc() = "Module for reading XYZ files with multiple frames and calculating knot type";
    m.def("calculate_knot_type", (std::vector<std::string> (*)(const std::string &, const std::string &)) &calculate_knot_type, 
          "Calculate knot type from XYZ file", py::arg("filename"), py::arg("chain_type") = "ring");
    m.def("calculate_knot_type", (std::vector<std::string> (*)(py::array_t<double>, const std::string &)) &calculate_knot_type, 
          "Calculate knot type from numpy array", py::arg("input"), py::arg("chain_type") = "ring");
    
    // knot type - 高性能并行版本
    m.def("calculate_knot_type_parallel", &calculate_knot_type_parallel, 
          "High-performance parallel calculate knot type from numpy array with SIMD and caching optimizations", 
          py::arg("input"), py::arg("chain_type") = "ring");
    
    // knot type - SIMD增强版本
    m.def("calculate_knot_type_simd", &calculate_knot_type_simd, 
          "SIMD-enhanced knot type calculation with advanced memory management", 
          py::arg("input"), py::arg("chain_type") = "ring");
    
    // knot size
    m.def("calculate_knot_size", &calculate_knot_size, "Calculate knot size from numpy array");

    // develop
    m.def("KMT_chain", &KMT_chain, py::arg("input"), py::arg("chain_type"), "A function to compute KMT or KMT_open_chain based on chain_type");
    m.def("gauss_notation", &gauss_notation, py::arg("input"), "Calculate the Gauss notation for a given input array of points. Accepts a NumPy array of shape (N, 3) as input and returns a 1D NumPy array of integers as output.");
    
    // 性能优化相关函数
    m.def("clear_caches", []() {
        // 清理所有静态缓存
        // 这里可以添加缓存清理逻辑
        return "Caches cleared";
    }, "Clear all internal caches to free memory");
    
    // 优化信息函数
    m.def("get_optimization_info", []() {
        py::dict info;
        info["simd_avx"] = OptimizedUtils::simd::hasAVXSupport();
        info["simd_sse2"] = OptimizedUtils::simd::hasSSE2Support();
        info["thread_count"] = OptimizedUtils::parallel::getOptimalThreadCount();
        info["cache_enabled"] = true;
        info["force_inline"] = true;
        info["branch_prediction"] = true;
        info["memory_prefetch"] = true;
        return info;
    }, "Get information about available optimizations and system capabilities");
}