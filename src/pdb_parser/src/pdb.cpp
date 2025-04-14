#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <limits>

namespace py = pybind11;

// 定义 Atom 结构体来存储原子信息
struct Atom {
    int serial;
    std::string name;
    std::string resName;
    std::string chainID;
    int resSeq; 
    double x, y, z;
    double occupancy;
    double tempFactor;
    std::string element;
    std::string charge;
    int frameIdx; // 新增：帧索引
};

// PDB 解析器类
class PDBParser {
public:
    PDBParser() = default;

    // 从文件加载 PDB 轨迹数据
    bool load(const std::string& filename) {
        std::ifstream file(filename);
        if (!file.is_open()) {
            std::cerr << "Can not open file: " << filename << std::endl;
            return false;
        }

        atoms.clear();
        frames = 0;
        
        std::string line;
        int currentFrame = 0;
        bool inNewFrame = true;
        
        while (std::getline(file, line)) {
            // 处理MODEL行，表示新帧开始
            if (line.substr(0, 5) == "MODEL") {
                inNewFrame = true;
                try {
                    // 提取MODEL后的所有数字内容，支持任意位数的编号
                    size_t startPos = line.find_first_not_of(" ", 5);
                    if (startPos != std::string::npos) {
                        std::string modelNumStr = line.substr(startPos);
                        // 清除尾部可能的空格
                        modelNumStr = trim(modelNumStr);
                        // 转换为整数并减1（因为MODEL编号通常从1开始）
                        currentFrame = std::stoi(modelNumStr) - 1;
                        if (currentFrame >= frames) frames = currentFrame + 1;
                    } else {
                        currentFrame = frames++;
                    }
                } catch (const std::exception& e) {
                    currentFrame = frames++;
                }
                continue;
            }
            
            // 处理ENDMDL行，表示当前帧结束
            if (line.substr(0, 6) == "ENDMDL") {
                inNewFrame = false;
                continue;
            }
            
            // 处理原子行
            if (line.substr(0, 4) == "ATOM" || line.substr(0, 6) == "HETATM") {
                // 确保没有MODEL标记的情况下也能正确处理单帧文件
                if (frames == 0) frames = 1;
                parseAtomLine(line, currentFrame);
            }
        }
        
        file.close();
        return true;
    }

    // 获取所有原子
    const std::vector<Atom>& getAtoms() const {
        return atoms;
    }
    
    // 获取特定帧的所有原子
    std::vector<Atom> getAtomsInFrame(int frameIdx) const {
        std::vector<Atom> frameAtoms;
        for (const auto& atom : atoms) {
            if (atom.frameIdx == frameIdx) {
                frameAtoms.push_back(atom);
            }
        }
        return frameAtoms;
    }

    // 获取原子坐标作为 numpy 数组
    py::array_t<double> getCoordinates() const {
        // 使用另一种构造方法
        std::vector<size_t> shape = {atoms.size(), 3};
        py::array_t<double> coords(shape);
        auto r = coords.mutable_unchecked<2>();
        
        for (size_t i = 0; i < atoms.size(); ++i) {
            r(i, 0) = atoms[i].x;
            r(i, 1) = atoms[i].y;
            r(i, 2) = atoms[i].z;
        }
        
        return coords;
    }

    // 获取特定帧的原子坐标
    py::array_t<double> getFrameCoordinates(int frameIdx) const {
        std::vector<Atom> frameAtoms = getAtomsInFrame(frameIdx);
        std::vector<size_t> shape = {frameAtoms.size(), 3};
        py::array_t<double> coords(shape);
        auto r = coords.mutable_unchecked<2>();
        
        for (size_t i = 0; i < frameAtoms.size(); ++i) {
            r(i, 0) = frameAtoms[i].x;
            r(i, 1) = frameAtoms[i].y;
            r(i, 2) = frameAtoms[i].z;
        }
        
        return coords;
    }
    
    // 获取所有CA原子坐标
    py::array_t<double> getCACoordinates() const {
        std::vector<Atom> caAtoms;
        
        for (const auto& atom : atoms) {
            // "CA  " 或 " CA " - 注意PDB中原子名称字段为4个字符，CA通常有空格填充
            std::string trimmedName = trim(atom.name);
            if (trimmedName == "CA") {
                caAtoms.push_back(atom);
            }
        }
        
        std::vector<size_t> shape = {caAtoms.size(), 3};
        py::array_t<double> coords(shape);
        auto r = coords.mutable_unchecked<2>();
        
        for (size_t i = 0; i < caAtoms.size(); ++i) {
            r(i, 0) = caAtoms[i].x;
            r(i, 1) = caAtoms[i].y;
            r(i, 2) = caAtoms[i].z;
        }
        
        return coords;
    }
    
    // 获取特定帧的CA原子坐标
    py::array_t<double> getFrameCACoordinates(int frameIdx) const {
        std::vector<Atom> caAtoms;
        
        for (const auto& atom : atoms) {
            std::string trimmedName = trim(atom.name);
            if (atom.frameIdx == frameIdx && trimmedName == "CA") {
                caAtoms.push_back(atom);
            }
        }
        
        std::vector<size_t> shape = {caAtoms.size(), 3};
        py::array_t<double> coords(shape);
        auto r = coords.mutable_unchecked<2>();
        
        for (size_t i = 0; i < caAtoms.size(); ++i) {
            r(i, 0) = caAtoms[i].x;
            r(i, 1) = caAtoms[i].y;
            r(i, 2) = caAtoms[i].z;
        }
        
        return coords;
    }
    
    // 获取帧数
    int getFrameCount() const {
        return frames;
    }

    py::array_t<double> getAllFrameCoordinates() const {
        if (frames == 0) return py::array_t<double>();
        
        // 先找出每一帧的原子数量
        std::vector<size_t> atomsPerFrame(frames, 0);
        for (const auto& atom : atoms) {
            if (atom.frameIdx >= 0 && atom.frameIdx < frames) {
                atomsPerFrame[atom.frameIdx]++;
            }
        }
        
        // 找出最大的原子数量，用于统一数组维度
        size_t maxAtoms = *std::max_element(atomsPerFrame.begin(), atomsPerFrame.end());
        
        // 创建3D数组 [frames, max_atoms, 3]
        std::vector<size_t> shape = {static_cast<size_t>(frames), maxAtoms, 3};
        py::array_t<double> allCoords(shape);
        auto r = allCoords.mutable_unchecked<3>();
        
        // 初始化所有值为NaN (表示此帧中没有该原子)
        for (size_t f = 0; f < frames; ++f) {
            for (size_t a = 0; a < maxAtoms; ++a) {
                for (size_t c = 0; c < 3; ++c) {
                    r(f, a, c) = std::numeric_limits<double>::quiet_NaN();
                }
            }
        }
        
        // 填充实际数据
        std::vector<size_t> atomCount(frames, 0);
        for (const auto& atom : atoms) {
            if (atom.frameIdx >= 0 && atom.frameIdx < frames) {
                size_t idx = atomCount[atom.frameIdx]++;
                if (idx < maxAtoms) {
                    r(atom.frameIdx, idx, 0) = atom.x;
                    r(atom.frameIdx, idx, 1) = atom.y;
                    r(atom.frameIdx, idx, 2) = atom.z;
                }
            }
        }
        
        return allCoords;
    }
    
    // 获取所有帧的CA原子坐标为一个3D数组
    py::array_t<double> getAllFrameCACoordinates() const {
        if (frames == 0) return py::array_t<double>();
        
        // 先计算每帧的CA原子数
        std::vector<std::vector<Atom>> frameCAAtoms(frames);
        for (const auto& atom : atoms) {
            if (atom.frameIdx >= 0 && atom.frameIdx < frames) {
                std::string trimmedName = trim(atom.name);
                if (trimmedName == "CA") {
                    frameCAAtoms[atom.frameIdx].push_back(atom);
                }
            }
        }
        
        // 找出最大CA原子数
        size_t maxCAAtoms = 0;
        for (const auto& frameAtoms : frameCAAtoms) {
            maxCAAtoms = std::max(maxCAAtoms, frameAtoms.size());
        }
        
        // 创建3D数组 [frames, max_ca_atoms, 3]
        std::vector<size_t> shape = {static_cast<size_t>(frames), maxCAAtoms, 3};
        py::array_t<double> allCACoords(shape);
        auto r = allCACoords.mutable_unchecked<3>();
        
        // 初始化所有值为NaN
        for (size_t f = 0; f < frames; ++f) {
            for (size_t a = 0; a < maxCAAtoms; ++a) {
                for (size_t c = 0; c < 3; ++c) {
                    r(f, a, c) = std::numeric_limits<double>::quiet_NaN();
                }
            }
        }
        
        // 填充实际数据
        for (size_t f = 0; f < frames; ++f) {
            const auto& frameAtoms = frameCAAtoms[f];
            for (size_t a = 0; a < frameAtoms.size(); ++a) {
                r(f, a, 0) = frameAtoms[a].x;
                r(f, a, 1) = frameAtoms[a].y;
                r(f, a, 2) = frameAtoms[a].z;
            }
        }
        
        return allCACoords;
    }

private:
    std::vector<Atom> atoms;
    int frames = 0;
    
    // 解析 ATOM/HETATM 行
    void parseAtomLine(const std::string& line, int frameIdx) {
        // PDB 文件格式解析 (根据 PDB 标准格式)
        if (line.size() < 54) return; // 最小需要坐标信息
        
        Atom atom;
        try {
            atom.serial = std::stoi(line.substr(6, 5));
            atom.name = line.substr(12, 4);
            atom.resName = line.substr(17, 3);
            atom.chainID = line.substr(21, 1);
            atom.resSeq = std::stoi(line.substr(22, 4));
            atom.x = std::stod(line.substr(30, 8));
            atom.y = std::stod(line.substr(38, 8));
            atom.z = std::stod(line.substr(46, 8));
            
            // 可选字段
            if (line.size() >= 60) atom.occupancy = std::stod(line.substr(54, 6));
            if (line.size() >= 66) atom.tempFactor = std::stod(line.substr(60, 6));
            if (line.size() >= 78) atom.element = line.substr(76, 2);
            if (line.size() >= 80) atom.charge = line.substr(78, 2);
            
            // 设置帧索引
            atom.frameIdx = frameIdx;
            
            atoms.push_back(atom);
        }
        catch (const std::exception& e) {
            std::cerr << "Parse error: " << e.what() << " in line: " << line << std::endl;
        }
    }
    
    // 辅助函数：去除字符串首尾空格
    static std::string trim(const std::string& str) {
        auto start = std::find_if_not(str.begin(), str.end(), [](char c) {
            return std::isspace(c);
        });
        
        auto end = std::find_if_not(str.rbegin(), str.rend(), [](char c) {
            return std::isspace(c);
        }).base();
        
        return (start < end) ? std::string(start, end) : std::string();
    }
};

// pybind11 模块定义
PYBIND11_MODULE(pdb_parser, m) {
    m.doc() = "PDB file parser for Python";
    
    py::class_<Atom>(m, "Atom")
        .def_readonly("serial", &Atom::serial)
        .def_readonly("name", &Atom::name)
        .def_readonly("resName", &Atom::resName)
        .def_readonly("chainID", &Atom::chainID)
        .def_readonly("resSeq", &Atom::resSeq)
        .def_readonly("x", &Atom::x)
        .def_readonly("y", &Atom::y)
        .def_readonly("z", &Atom::z)
        .def_readonly("occupancy", &Atom::occupancy)
        .def_readonly("tempFactor", &Atom::tempFactor)
        .def_readonly("element", &Atom::element)
        .def_readonly("charge", &Atom::charge)
        .def_readonly("frameIdx", &Atom::frameIdx)
        .def("__repr__", [](const Atom& a) {
            return "<Atom " + std::to_string(a.serial) + ": " + a.name + " (Frame: " + 
                   std::to_string(a.frameIdx) + ", coordinates: " + 
                   std::to_string(a.x) + ", " + std::to_string(a.y) + ", " + 
                   std::to_string(a.z) + ")>";
        });
    
    py::class_<PDBParser>(m, "PDBParser")
        .def(py::init<>())
        .def("load", &PDBParser::load, "从文件加载PDB轨迹数据")
        .def("get_atoms", &PDBParser::getAtoms, "获取所有原子")
        .def("get_atoms_in_frame", &PDBParser::getAtomsInFrame, "获取特定帧的所有原子")
        .def("get_coordinates", &PDBParser::getCoordinates, "获取所有原子坐标作为numpy数组")
        .def("get_frame_coordinates", &PDBParser::getFrameCoordinates, "获取特定帧的原子坐标")
        .def("get_ca_coordinates", &PDBParser::getCACoordinates, "获取所有CA原子坐标")
        .def("get_frame_ca_coordinates", &PDBParser::getFrameCACoordinates, "获取特定帧的CA原子坐标")
        .def("get_frame_count", &PDBParser::getFrameCount, "获取PDB文件中的帧数")
        .def("get_all_frame_coordinates", &PDBParser::getAllFrameCoordinates, "获取所有帧的原子坐标")
        .def("get_all_frame_ca_coordinates", &PDBParser::getAllFrameCACoordinates, "获取所有帧的CA原子坐标");
}
