# 性能测试使用指南

本项目提供了两个性能测试脚本来评估优化效果：

## 📊 测试文件说明

### 1. `benchmark_performance.py` - 基础性能测试
- **用途**: 测试核心函数的基本性能指标
- **特点**: 多次运行取平均值，包含预热阶段
- **适用场景**: 日常性能监控和回归测试

### 2. `test_optimized_performance.py` - 高级优化测试
- **用途**: 测试缓存、内存池等优化特性
- **特点**: 专门测试缓存命中率、并发性能等
- **适用场景**: 验证优化效果和深度性能分析

## 🚀 使用方法

### 步骤1: 编译C++模块

```bash
# 进入C++模块目录
cd src/cpp_module

# 编译模块
python setup.py build_ext --inplace

# 或者使用cmake编译
mkdir -p build && cd build
cmake ..
make -j$(nproc)
```

### 步骤2: 确保测试数据存在

确保以下测试文件存在于 `test/` 目录：
- `traj_knot31_L300_close_short.txt`
- `traj_knot31_L300_open_short.txt`

### 步骤3: 运行基础性能测试

```bash
# 在项目根目录运行
cd /home/yongjian/Yongjian_data/github/pythonknot
python test/benchmark_performance.py
```

**预期输出示例:**
```
==========================================================
PERFORMANCE BENCHMARK TESTS
==========================================================
✓ Successfully imported alexander_poly module
Loading test data...
Loaded data with shape: (50, 300, 3)

--- Testing with single frame (1 frames, 300 atoms) ---
Benchmarking calculate_knot_type (ring):
  Mean time: 0.0234s ± 0.0012s
  Range: 0.0221s - 0.0248s
  Results: 1 knot types

Benchmarking calculate_knot_type (open):
  Mean time: 0.0198s ± 0.0008s
  Range: 0.0189s - 0.0207s
  Results: 1 knot types
```

### 步骤4: 运行优化测试

```bash
python test/test_optimized_performance.py
```

**预期输出示例:**
```
开始优化性能测试...
初始内存使用: 45.23 MB
模块验证成功，测试结果: ['0_1']

=== 缓存性能测试 ===
测试数据形状: (50, 300, 3)
重复数据形状: (500, 300, 3)
第一次运行（冷启动）:
  时间: 2.3456s
  内存增长: 12.34 MB
第二次运行（热启动，应使用缓存）:
  时间: 0.1234s
  内存增长: 0.56 MB
缓存加速比: 19.01x
```

## 📈 性能指标解读

### 时间指标
- **平均时间**: 多次运行的平均执行时间
- **标准差**: 时间波动程度，越小越稳定
- **范围**: 最快和最慢的执行时间

### 内存指标
- **内存增长**: 函数执行前后的内存差异
- **内存峰值**: 执行过程中的最大内存使用
- **缓存效率**: 第二次运行相比第一次的提升

### 缓存指标
- **缓存命中率**: 相同数据被缓存命中的比例
- **加速比**: 缓存生效后的性能提升倍数
- **内存节省**: 缓存避免的重复计算开销

## 🔧 故障排除

### 常见问题1: 模块导入失败
```
ImportError: No module named 'alexander_poly'
```
**解决方案:**
```bash
# 检查模块是否编译成功
ls src/cpp_module/build/
# 应该看到 alexander_poly.cpython-*.so 文件

# 重新编译
cd src/cpp_module
python setup.py build_ext --inplace --force
```

### 常见问题2: 测试数据缺失
```
FileNotFoundError: [Errno 2] No such file or directory: 'test/traj_knot31_L300_close_short.txt'
```
**解决方案:**
```bash
# 检查测试数据
ls test/*.txt
# 确保存在必要的测试文件
```

### 常见问题3: C++17兼容性错误
```
error: 'if constexpr' only available with -std=c++17
```
**解决方案:**
- 当前代码已优化为C++11兼容
- 确保使用最新版本的优化代码

## 📊 性能基线参考

基于300原子、50帧的测试数据：

| 函数 | 单帧时间 | 内存使用 | 缓存加速 |
|------|----------|----------|----------|
| calculate_knot_type | ~0.02s | ~2MB | 10-20x |
| calculate_knot_size | ~0.05s | ~4MB | 5-10x |
| read_xyz | ~0.01s | ~5MB | N/A |
| KMT_chain | ~0.03s | ~3MB | N/A |

## 🎯 优化验证清单

运行测试后，检查以下指标：

- [ ] **缓存加速比 > 5x**: 重复数据处理显著加速
- [ ] **内存增长 < 10MB**: 单次操作内存控制合理
- [ ] **时间标准差 < 10%**: 性能稳定性良好
- [ ] **并发测试通过**: 多线程安全性确认
- [ ] **大数据集可处理**: 1000+帧数据不崩溃

如果某项指标不达标，请检查编译配置和代码版本。

## 🔍 高级分析

### 使用profiler深度分析
```bash
# 安装profiler
pip install line_profiler memory_profiler

# 行级性能分析
kernprof -l -v test/benchmark_performance.py

# 内存使用分析
mprof run test/test_optimized_performance.py
mprof plot
```

### 自定义测试场景
```python
# 在benchmark_performance.py中添加自定义测试
def custom_benchmark():
    # 你的测试数据
    custom_data = np.random.rand(100, 500, 3)
    
    # 测试特定函数
    stats = benchmark_function(alexander_poly.calculate_knot_type, 
                              custom_data, "ring")
    print(f"Custom test: {stats['mean']:.4f}s")
```

通过这些测试，你可以全面评估优化效果并确保代码质量！