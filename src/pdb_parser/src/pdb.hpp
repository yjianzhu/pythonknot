#pragma once

#include <string>
#include <vector>
#include <pybind11/numpy.h>

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
    PDBParser();
    
    // 从文件加载 PDB 轨迹数据
    bool load(const std::string& filename);
    
    // 获取所有原子
    const std::vector<Atom>& getAtoms() const;
    
    // 获取特定帧的所有原子
    std::vector<Atom> getAtomsInFrame(int frameIdx) const;
    
    // 获取原子坐标作为 numpy 数组
    py::array_t<double> getCoordinates() const;
    
    // 获取特定帧的原子坐标
    py::array_t<double> getFrameCoordinates(int frameIdx) const;
    
    // 获取所有CA原子坐标
    py::array_t<double> getCACoordinates() const;
    
    // 获取特定帧的CA原子坐标
    py::array_t<double> getFrameCACoordinates(int frameIdx) const;

    py::array_t<double> getAllFrameCoordinates() const;

    py::array_t<double> getAllFrameCACoordinates() const;
    
    // 获取帧数
    int getFrameCount() const;

private:
    std::vector<Atom> atoms;
    int frames = 0;
    
    // 解析 ATOM/HETATM 行
    void parseAtomLine(const std::string& line, int frameIdx);
    
    // 辅助函数：去除字符串首尾空格
    static std::string trim(const std::string& str);
};