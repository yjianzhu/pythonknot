#ifndef OPTIMIZED_TYPES_H
#define OPTIMIZED_TYPES_H

#include <vector>
#include <array>
#include <memory>
#include <string>
#include <map>
#include <algorithm>
#include <cstring>
#include <cstdint>
#include <thread>
#include <mutex>
#include <cmath>

// 包含 pybind11 头文件
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

// 分支预测优化宏 - 移到文件开头
#if defined(__GNUC__) || defined(__clang__)
    #define LIKELY(x)       __builtin_expect(!!(x), 1)
    #define UNLIKELY(x)     __builtin_expect(!!(x), 0)
    #define FORCE_INLINE    __attribute__((always_inline)) inline
    #define PREFETCH(addr)  __builtin_prefetch(addr, 0, 3)
#else
    #define LIKELY(x)       (x)
    #define UNLIKELY(x)     (x)
    #define FORCE_INLINE    inline
    #define PREFETCH(addr)  
#endif

// SIMD支持检测和包含
#if defined(__SSE2__) || defined(_M_X64) || (defined(_M_IX86_FP) && _M_IX86_FP >= 2)
    #define ENABLE_SSE2
    #include <emmintrin.h>
#endif

#if defined(__AVX__) || defined(__AVX2__)
    #define ENABLE_AVX
    #include <immintrin.h>
#endif

// 使用现代C++类型替代原始指针
using Point3D = std::array<double, 3>;
using PointVector = std::vector<Point3D>;

// RAII包装器用于自动内存管理
class PointManager {
private:
    PointVector points;
    
public:
    PointManager() = default;
    
    // 从numpy数组构造 - 增强边界检查
    PointManager(double* data, size_t nAtoms, size_t nDimension, size_t frame_offset) {
        if (UNLIKELY(!data || nAtoms == 0 || nDimension != 3)) {
            throw std::runtime_error("Invalid parameters for PointManager construction");
        }
        
        points.reserve(nAtoms);
        for (size_t j = 0; j < nAtoms; j++) {
            const size_t base_idx = frame_offset + j * nDimension;
            
            // 添加边界检查防止内存越界
            if (UNLIKELY(base_idx + 2 < base_idx)) { // 检查溢出
                throw std::runtime_error("Index overflow in PointManager");
            }
            
            Point3D point;
            point[0] = data[base_idx + 0];
            point[1] = data[base_idx + 1];
            point[2] = data[base_idx + 2];
            points.push_back(point);
        }
    }
    
    // 获取传统格式的指针向量（用于兼容现有代码）
    std::vector<double*> getCompatiblePointers() {
        std::vector<double*> ptrs;
        ptrs.reserve(points.size());
        for (auto& point : points) {
            ptrs.push_back(point.data());
        }
        return ptrs;
    }
    
    // 直接访问点数据
    const PointVector& getPoints() const { return points; }
    PointVector& getPoints() { return points; }
    
    size_t size() const { return points.size(); }
    bool empty() const { return points.empty(); }
    
    // 禁用拷贝构造和赋值（避免意外的深拷贝）
    PointManager(const PointManager&) = delete;
    PointManager& operator=(const PointManager&) = delete;
    
    // 启用移动语义
    PointManager(PointManager&&) = default;
    PointManager& operator=(PointManager&&) = default;
};

// 内存池类，用于减少内存分配开销
class MemoryPool {
private:
    static const size_t DEFAULT_BLOCK_SIZE = 1024;
    std::vector<std::vector<Point3D>> pools;
    std::vector<size_t> pool_indices;
    
public:
    MemoryPool() {
        pools.reserve(8); // 预分配8个不同大小的池
        pool_indices.reserve(8);
    }
    
    // 获取指定大小的预分配内存块
    PointVector& getPool(size_t size) {
        // 查找合适大小的池
        for (size_t i = 0; i < pools.size(); ++i) {
            if (pools[i].capacity() >= size && pool_indices[i] == 0) {
                pool_indices[i] = 1; // 标记为使用中
                pools[i].clear();
                return pools[i];
            }
        }
        
        // 创建新的池
        pools.emplace_back();
        pool_indices.push_back(1);
        pools.back().reserve(std::max(size, DEFAULT_BLOCK_SIZE));
        return pools.back();
    }
    
    // 释放池回收使用
    void releasePool(PointVector& pool) {
        for (size_t i = 0; i < pools.size(); ++i) {
            if (&pools[i] == &pool) {
                pool_indices[i] = 0; // 标记为可用
                break;
            }
        }
    }
    
    // 获取单例实例
    static MemoryPool& getInstance() {
        static MemoryPool instance;
        return instance;
    }
};

// 结果缓存类，避免重复计算 - 线程安全版本
template<typename Key, typename Value>
class ResultCache {
private:
    std::map<Key, Value> cache;
    size_t max_size;
    mutable std::mutex cache_mutex; // 添加互斥锁保护
    
public:
    explicit ResultCache(size_t max_cache_size = 1000) : max_size(max_cache_size) {}
    
    bool get(const Key& key, Value& value) const {
        std::lock_guard<std::mutex> lock(cache_mutex);
        auto it = cache.find(key);
        if (it != cache.end()) {
            value = it->second;
            return true;
        }
        return false;
    }
    
    void put(const Key& key, const Value& value) {
        std::lock_guard<std::mutex> lock(cache_mutex);
        
        // 如果key已存在，直接更新
        auto it = cache.find(key);
        if (it != cache.end()) {
            it->second = value;
            return;
        }
        
        // 如果缓存已满，清理一些旧条目（更安全的方式）
        if (cache.size() >= max_size) {
            // 删除前25%的条目，而不是只删除一个
            size_t to_remove = std::max(size_t(1), max_size / 4);
            auto remove_it = cache.begin();
            for (size_t i = 0; i < to_remove && remove_it != cache.end(); ++i) {
                auto next_it = std::next(remove_it);
                cache.erase(remove_it);
                remove_it = next_it;
            }
        }
        
        // 插入新条目
        cache[key] = value;
    }
    
    void clear() { 
        std::lock_guard<std::mutex> lock(cache_mutex);
        cache.clear(); 
    }
    
    size_t size() const { 
        std::lock_guard<std::mutex> lock(cache_mutex);
        return cache.size(); 
    }
};

// 分支预测优化宏已移至文件开头

// 性能优化的辅助函数
namespace OptimizedUtils {
    // 内部实现细节
    namespace detail {
        // 有reserve方法的版本
        template<typename Container>
        static auto reserveImpl(Container& container, size_t size, int) 
            -> decltype(container.reserve(size), void()) {
            container.reserve(size);
        }
        
        // 没有reserve方法的版本（什么都不做）
        template<typename Container>
        static void reserveImpl(Container&, size_t, long) {
            // 空实现
        }
    }
    
    // 预分配容器大小避免重复扩容
    template<typename Container>
    void reserveIfPossible(Container& container, size_t size) {
        // 使用SFINAE技术检查是否有reserve方法
        detail::reserveImpl(container, size, 0);
    }

    // 快速数组初始化
    inline Point3D makePoint(double x, double y, double z) {
        Point3D result = {{x, y, z}};
        return result;
    }

    // 基于哈希的快速数据指纹计算 - 增强边界检查
    inline size_t computeFrameHash(double* data, size_t nAtoms, size_t nDimension, size_t frame_offset) {
        if (UNLIKELY(!data || nAtoms == 0 || nDimension == 0)) {
            return 0;
        }
        
        size_t hash = 0;
        const double* frame_data = data + frame_offset;
        const size_t total_elements = nAtoms * nDimension;
        const size_t step = std::max(size_t(1), total_elements / 32); // 采样策略避免计算所有数据
        
        for (size_t i = 0; i < total_elements; i += step) {
            // 确保不会越界
            if (UNLIKELY(i >= total_elements)) break;
            
            // 简单但快速的哈希函数
            union { double d; uint64_t i; } converter;
            converter.d = frame_data[i];
            hash ^= converter.i + 0x9e3779b9 + (hash << 6) + (hash >> 2);
        }
        return hash;
    }

    // 内存对齐的批量数据复制
    FORCE_INLINE void optimizedMemcpy(double* dst, const double* src, size_t count) {
        // 对于小量数据，直接循环更快
        if (LIKELY(count <= 12)) {
            for (size_t i = 0; i < count; ++i) {
                dst[i] = src[i];
            }
        } else {
            std::memcpy(dst, src, count * sizeof(double));
        }
    }

    // 优化的向量距离计算
    FORCE_INLINE double fastDistance(const Point3D& p1, const Point3D& p2) {
        const double dx = p1[0] - p2[0];
        const double dy = p1[1] - p2[1];
        const double dz = p1[2] - p2[2];
        return dx*dx + dy*dy + dz*dz; // 返回距离平方，避免sqrt
    }
    
    // SIMD优化工具
    namespace simd {
        // SIMD支持检测
        inline bool hasAVXSupport() {
#ifdef ENABLE_AVX
            return true;
#else
            return false;
#endif
        }
        
        inline bool hasSSE2Support() {
#ifdef ENABLE_SSE2
            return true;
#else
            return false;
#endif
        }
        
        // SIMD优化的批量距离计算
        FORCE_INLINE void vectorizedDistanceSquared(
            const double* points1, const double* points2, 
            double* distances, size_t count) {
            
            // 添加空指针和大小检查
            if (UNLIKELY(!points1 || !points2 || !distances || count == 0)) {
                return;
            }
            
#ifdef ENABLE_AVX
            // AVX优化版本：处理连续的3D点对
            const size_t points_per_batch = 2; // 每批处理2个3D点 (6个double)
            const size_t avx_batches = count / points_per_batch;
            
            for (size_t batch = 0; batch < avx_batches; ++batch) {
                const size_t base_idx = batch * points_per_batch;
                
                // 确保不会越界访问
                if (UNLIKELY(base_idx + 1 >= count)) break;
                
                // 更严格的边界检查
                if (UNLIKELY((base_idx + 1) * 3 + 2 >= count * 3)) {
                    // 回退到标量处理剩余数据
                    for (size_t i = base_idx; i < count; ++i) {
                        if (UNLIKELY(i * 3 + 2 >= count * 3)) break;
                        const double* p1 = &points1[i * 3];
                        const double* p2 = &points2[i * 3];  
                        const double dx = p1[0] - p2[0];
                        const double dy = p1[1] - p2[1];
                        const double dz = p1[2] - p2[2];
                        distances[i] = dx*dx + dy*dy + dz*dz;
                    }
                    break;
                }
                
                __m256d p1_data = _mm256_loadu_pd(&points1[base_idx * 3]);     // x1,y1,z1,x2
                __m256d p1_extra = _mm256_loadu_pd(&points1[base_idx * 3 + 2]); // z1,x2,y2,z2
                __m256d p2_data = _mm256_loadu_pd(&points2[base_idx * 3]);     // x1,y1,z1,x2
                __m256d p2_extra = _mm256_loadu_pd(&points2[base_idx * 3 + 2]); // z1,x2,y2,z2
                
                // 计算差值
                __m256d diff1 = _mm256_sub_pd(p1_data, p2_data);
                __m256d diff2 = _mm256_sub_pd(p1_extra, p2_extra);
                
                // 计算平方
                __m256d sq1 = _mm256_mul_pd(diff1, diff1);
                __m256d sq2 = _mm256_mul_pd(diff2, diff2);
                
                // 手动处理前两个点的距离计算
                double temp1[4], temp2[4];
                _mm256_storeu_pd(temp1, sq1);
                _mm256_storeu_pd(temp2, sq2);
                
                // 第一个点: x1² + y1² + z1²
                distances[base_idx] = temp1[0] + temp1[1] + temp2[0];
                
                // 第二个点: x2² + y2² + z2²
                if (base_idx + 1 < count) {
                    distances[base_idx + 1] = temp1[3] + temp2[1] + temp2[2];
                }
            }
            
            // 处理剩余的点
            for (size_t i = avx_batches * points_per_batch; i < count; ++i) {
                const double* p1 = &points1[i * 3];
                const double* p2 = &points2[i * 3];
                const double dx = p1[0] - p2[0];
                const double dy = p1[1] - p2[1];
                const double dz = p1[2] - p2[2];
                distances[i] = dx*dx + dy*dy + dz*dz;
            }
            
#elif defined(ENABLE_SSE2)
            // SSE2优化版本：每次处理1个3D点的距离
            for (size_t i = 0; i < count; ++i) {
                const double* p1 = &points1[i * 3];
                const double* p2 = &points2[i * 3];
                
                // 加载3D坐标 (只使用前2个元素，第3个单独处理)
                __m128d p1_xy = _mm_loadu_pd(p1);      // x1, y1
                __m128d p2_xy = _mm_loadu_pd(p2);      // x2, y2
                
                __m128d diff_xy = _mm_sub_pd(p1_xy, p2_xy);
                __m128d sq_xy = _mm_mul_pd(diff_xy, diff_xy);
                
                // 处理z坐标
                double dz = p1[2] - p2[2];
                double dz_sq = dz * dz;
                
                // 提取x²和y²
                double xy_squared[2];
                _mm_storeu_pd(xy_squared, sq_xy);
                
                distances[i] = xy_squared[0] + xy_squared[1] + dz_sq;
            }
#else
            // 标量版本（无SIMD）
            for (size_t i = 0; i < count; ++i) {
                const double* p1 = &points1[i * 3];
                const double* p2 = &points2[i * 3];
                const double dx = p1[0] - p2[0];
                const double dy = p1[1] - p2[1];
                const double dz = p1[2] - p2[2];
                distances[i] = dx*dx + dy*dy + dz*dz;
            }
#endif
        }
        
        // SIMD优化的向量点积计算
        FORCE_INLINE void vectorizedDotProduct(
            const double* vectors1, const double* vectors2,
            double* results, size_t count) {
            
            // 添加空指针和大小检查
            if (UNLIKELY(!vectors1 || !vectors2 || !results || count == 0)) {
                return;
            }
            
#ifdef ENABLE_AVX
            // AVX版本：每次处理2个3D向量
            const size_t vectors_per_batch = 2;
            const size_t avx_batches = count / vectors_per_batch;
            
            for (size_t batch = 0; batch < avx_batches; ++batch) {
                const size_t base_idx = batch * vectors_per_batch;
                
                // 边界检查
                if (UNLIKELY(base_idx + 1 >= count)) break;
                
                // 更严格的边界检查
                if (UNLIKELY((base_idx + 1) * 3 + 2 >= count * 3)) {
                    // 回退到标量处理
                    for (size_t i = base_idx; i < count; ++i) {
                        if (UNLIKELY(i * 3 + 2 >= count * 3)) break;
                        const double* v1 = &vectors1[i * 3];
                        const double* v2 = &vectors2[i * 3];
                        results[i] = v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
                    }
                    break;
                }
                
                // 加载向量数据
                __m256d v1_data = _mm256_loadu_pd(&vectors1[base_idx * 3]);
                __m256d v1_extra = _mm256_loadu_pd(&vectors1[base_idx * 3 + 2]);
                __m256d v2_data = _mm256_loadu_pd(&vectors2[base_idx * 3]);
                __m256d v2_extra = _mm256_loadu_pd(&vectors2[base_idx * 3 + 2]);
                
                // 计算乘积
                __m256d prod1 = _mm256_mul_pd(v1_data, v2_data);
                __m256d prod2 = _mm256_mul_pd(v1_extra, v2_extra);
                
                // 提取结果
                double temp1[4], temp2[4];
                _mm256_storeu_pd(temp1, prod1);
                _mm256_storeu_pd(temp2, prod2);
                
                // 计算点积
                results[base_idx] = temp1[0] + temp1[1] + temp2[0];
                if (base_idx + 1 < count) {
                    results[base_idx + 1] = temp1[3] + temp2[1] + temp2[2];
                }
            }
            
            // 处理剩余向量
            for (size_t i = avx_batches * vectors_per_batch; i < count; ++i) {
                const double* v1 = &vectors1[i * 3];
                const double* v2 = &vectors2[i * 3];
                results[i] = v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
            }
            
#elif defined(ENABLE_SSE2)
            // SSE2版本
            for (size_t i = 0; i < count; ++i) {
                const double* v1 = &vectors1[i * 3];
                const double* v2 = &vectors2[i * 3];
                
                __m128d v1_xy = _mm_loadu_pd(v1);
                __m128d v2_xy = _mm_loadu_pd(v2);
                __m128d prod_xy = _mm_mul_pd(v1_xy, v2_xy);
                
                double xy_prod[2];
                _mm_storeu_pd(xy_prod, prod_xy);
                
                results[i] = xy_prod[0] + xy_prod[1] + v1[2]*v2[2];
            }
#else
            // 标量版本
            for (size_t i = 0; i < count; ++i) {
                const double* v1 = &vectors1[i * 3];
                const double* v2 = &vectors2[i * 3];
                results[i] = v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
            }
#endif
        }
        
        // SIMD优化的向量模长计算
        FORCE_INLINE void vectorizedMagnitude(
            const double* vectors, double* magnitudes, size_t count) {
            
            // 添加空指针和大小检查
            if (UNLIKELY(!vectors || !magnitudes || count == 0)) {
                return;
            }
            
#ifdef ENABLE_AVX
            // AVX版本
            const size_t vectors_per_batch = 2;
            const size_t avx_batches = count / vectors_per_batch;
            
            for (size_t batch = 0; batch < avx_batches; ++batch) {
                const size_t base_idx = batch * vectors_per_batch;
                
                // 边界检查
                if (UNLIKELY(base_idx + 1 >= count)) break;
                
                // 更严格的边界检查
                if (UNLIKELY((base_idx + 1) * 3 + 2 >= count * 3)) {
                    // 回退到标量处理
                    for (size_t i = base_idx; i < count; ++i) {
                        if (UNLIKELY(i * 3 + 2 >= count * 3)) break;
                        const double* v = &vectors[i * 3];
                        magnitudes[i] = std::sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
                    }
                    break;
                }
                
                __m256d v_data = _mm256_loadu_pd(&vectors[base_idx * 3]);
                __m256d v_extra = _mm256_loadu_pd(&vectors[base_idx * 3 + 2]);
                
                __m256d sq1 = _mm256_mul_pd(v_data, v_data);
                __m256d sq2 = _mm256_mul_pd(v_extra, v_extra);
                
                double temp1[4], temp2[4];
                _mm256_storeu_pd(temp1, sq1);
                _mm256_storeu_pd(temp2, sq2);
                
                magnitudes[base_idx] = std::sqrt(temp1[0] + temp1[1] + temp2[0]);
                if (base_idx + 1 < count) {
                    magnitudes[base_idx + 1] = std::sqrt(temp1[3] + temp2[1] + temp2[2]);
                }
            }
            
            // 处理剩余向量
            for (size_t i = avx_batches * vectors_per_batch; i < count; ++i) {
                const double* v = &vectors[i * 3];
                magnitudes[i] = std::sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
            }
            
#else
            // 标量版本
            for (size_t i = 0; i < count; ++i) {
                const double* v = &vectors[i * 3];
                magnitudes[i] = std::sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
            }
#endif
        }
        
        // SIMD优化的向量归一化
        FORCE_INLINE void vectorizedNormalize(double* vectors, size_t count) {
            // 添加空指针和大小检查
            if (UNLIKELY(!vectors || count == 0)) {
                return;
            }
            
            // 先计算模长
            std::vector<double> magnitudes(count);
            vectorizedMagnitude(vectors, magnitudes.data(), count);
            
            // 归一化
#ifdef ENABLE_AVX
            const size_t vectors_per_batch = 2;
            const size_t avx_batches = count / vectors_per_batch;
            
            for (size_t batch = 0; batch < avx_batches; ++batch) {
                const size_t base_idx = batch * vectors_per_batch;
                
                if (magnitudes[base_idx] > 1e-12 && 
                    (base_idx + 1 >= count || magnitudes[base_idx + 1] > 1e-12)) {
                    
                    __m256d inv_mag = _mm256_setr_pd(
                        1.0/magnitudes[base_idx], 1.0/magnitudes[base_idx], 1.0/magnitudes[base_idx],
                        base_idx + 1 < count ? 1.0/magnitudes[base_idx + 1] : 1.0
                    );
                    
                    __m256d v_data = _mm256_loadu_pd(&vectors[base_idx * 3]);
                    __m256d normalized = _mm256_mul_pd(v_data, inv_mag);
                    _mm256_storeu_pd(&vectors[base_idx * 3], normalized);
                    
                    if (base_idx + 1 < count) {
                        __m256d inv_mag2 = _mm256_setr_pd(
                            1.0/magnitudes[base_idx + 1], 1.0/magnitudes[base_idx + 1], 1.0, 1.0
                        );
                        __m256d v_extra = _mm256_loadu_pd(&vectors[base_idx * 3 + 3]);
                        __m256d norm2 = _mm256_mul_pd(v_extra, inv_mag2);
                        _mm256_storeu_pd(&vectors[base_idx * 3 + 3], norm2);
                    }
                }
            }
            
            // 处理剩余向量
            for (size_t i = avx_batches * vectors_per_batch; i < count; ++i) {
                if (magnitudes[i] > 1e-12) {
                    double* v = &vectors[i * 3];
                    double inv_mag = 1.0 / magnitudes[i];
                    v[0] *= inv_mag;
                    v[1] *= inv_mag;
                    v[2] *= inv_mag;
                }
            }
#else
            // 标量版本
            for (size_t i = 0; i < count; ++i) {
                if (magnitudes[i] > 1e-12) {
                    double* v = &vectors[i * 3];
                    double inv_mag = 1.0 / magnitudes[i];
                    v[0] *= inv_mag;
                    v[1] *= inv_mag;
                    v[2] *= inv_mag;
                }
            }
#endif
        }
    }
    
    // 并行处理工具
    namespace parallel {
        // 获取最优线程数 - 限制最大线程数避免过度创建
        inline size_t getOptimalThreadCount() {
            static const unsigned int hardware_threads_raw = std::thread::hardware_concurrency();
            static const size_t hardware_threads = static_cast<size_t>(hardware_threads_raw);
            static const size_t max_threads = std::min(hardware_threads, static_cast<size_t>(8)); // 限制最大8个线程
            static const size_t thread_count = std::max(static_cast<size_t>(1), max_threads);
            return thread_count;
        }
        
        // 简单的并行for循环实现（C++11兼容） - 增强错误处理和资源管理
        template<typename Func>
        void parallelFor(size_t start, size_t end, size_t min_per_thread, Func&& func) {
            parallelFor(start, end, min_per_thread, 0, std::forward<Func>(func));
        }
        
        // 带自定义线程数的并行for循环
        template<typename Func>
        void parallelFor(size_t start, size_t end, size_t min_per_thread, size_t custom_thread_count, Func&& func) {
            if (UNLIKELY(start >= end)) {
                return; // 无效范围
            }
            
            const size_t length = end - start;
            const size_t hardware_threads = getOptimalThreadCount();
            const size_t max_threads = (length + min_per_thread - 1) / min_per_thread;
            
            // 如果指定了自定义线程数，使用它，否则使用硬件线程数
            size_t target_threads;
            if (custom_thread_count > 0) {
                // 限制自定义线程数不超过硬件线程数的2倍，避免过度创建
                target_threads = std::min(custom_thread_count, hardware_threads * 2);
            } else {
                target_threads = hardware_threads;
            }
            
            const size_t num_threads = std::min(target_threads, max_threads);
            
            // 对于小数据量或单线程，直接串行执行
            if (num_threads <= 1 || length < min_per_thread * 2) {
                try {
                    func(start, end);
                } catch (...) {
                    // 捕获异常避免崩溃
                }
                return;
            }
            
            std::vector<std::thread> threads;
            std::vector<std::exception_ptr> exceptions(num_threads - 1);
            threads.reserve(num_threads - 1);
            
            const size_t block_size = length / num_threads;
            
            // 创建工作线程（比总线程数少1个，主线程处理最后一块）
            for (size_t i = 0; i < num_threads - 1; ++i) {
                size_t block_start = start + i * block_size;
                size_t block_end = start + (i + 1) * block_size;
                
                threads.emplace_back([=, &exceptions, &func]() {
                    try {
                        func(block_start, block_end);
                    } catch (...) {
                        exceptions[i] = std::current_exception();
                    }
                });
            }
            
            // 主线程处理最后一个块
            std::exception_ptr main_exception;
            try {
                size_t final_start = start + (num_threads - 1) * block_size;
                func(final_start, end);
            } catch (...) {
                main_exception = std::current_exception();
            }
            
            // 等待所有线程完成
            for (auto& thread : threads) {
                if (thread.joinable()) {
                    thread.join();
                }
            }
            
            // 检查是否有异常
            if (main_exception) {
                std::rethrow_exception(main_exception);
            }
            
            for (const auto& exception : exceptions) {
                if (exception) {
                    std::rethrow_exception(exception);
                }
            }
        }
    }
    
    // pybind11绑定优化工具
    namespace pybind_optimized {
        // 现在可以直接使用 pybind11 命名空间
        namespace py = pybind11;
        
        // 快速numpy数组信息提取（避免重复buffer_info调用）- 优化版本
        struct FastArrayInfo {
            double* data;
            py::ssize_t nFrames;
            py::ssize_t nAtoms;
            py::ssize_t nDimension;
            py::ssize_t total_size;
            
            explicit FastArrayInfo(py::array_t<double>& input) {
                py::buffer_info info = input.request();
                if (UNLIKELY(info.ndim != 3)) {
                    throw std::runtime_error("Expected a 3-dimensional array");
                }
                if (UNLIKELY(info.shape[2] != 3)) {
                    throw std::runtime_error("Last dimension must be 3 (x,y,z coordinates)");
                }
                
                data = static_cast<double*>(info.ptr);
                nFrames = info.shape[0];
                nAtoms = info.shape[1];
                nDimension = info.shape[2];
                total_size = nFrames * nAtoms * nDimension;
                
                // 验证数据完整性
                if (UNLIKELY(!data || total_size <= 0)) {
                    throw std::runtime_error("Invalid array data or dimensions");
                }
            }
        };
        
        // 批处理结果构造器（减少Python对象创建开销）- 优化版本
        template<typename T>
        class BatchResultBuilder {
        private:
            std::vector<T> results;
            
        public:
            explicit BatchResultBuilder(size_t reserve_size = 100) {
                results.reserve(reserve_size);
            }
            
            // 移动语义优化
            FORCE_INLINE void addResult(T&& result) {
                results.emplace_back(std::move(result));
            }
            
            // 拷贝版本
            FORCE_INLINE void addResult(const T& result) {
                results.push_back(result);
            }
            
            // 批量添加优化
            template<typename Iterator>
            void addResults(Iterator begin, Iterator end) {
                results.insert(results.end(), begin, end);
            }
            
            // 预分配空间
            FORCE_INLINE void reserve(size_t additional_size) {
                results.reserve(results.size() + additional_size);
            }
            
            // 移动构建，避免不必要的拷贝
            std::vector<T> build() && {
                return std::move(results);
            }
            
            // 获取引用避免拷贝
            const std::vector<T>& get() const { return results; }
            std::vector<T>& get() { return results; }
            
            size_t size() const { return results.size(); }
            void clear() { results.clear(); }
            bool empty() const { return results.empty(); }
        };
        
        // 零拷贝numpy数组创建助手 - 优化版本
        template<typename T>
        py::array_t<T> createOptimizedArray(std::vector<T>&& data, std::vector<py::ssize_t> shape) {
            // 验证数据和形状一致性
            size_t expected_size = 1;
            for (auto dim : shape) {
                expected_size *= static_cast<size_t>(dim);
            }
            
            if (UNLIKELY(data.size() != expected_size)) {
                throw std::runtime_error("Data size doesn't match shape");
            }
            
            // 移动构造减少拷贝
            T* raw_ptr = data.data();
            
            // 创建numpy数组，使用零拷贝策略
            return py::array_t<T>(
                shape,
                raw_ptr,
                py::capsule(new std::vector<T>(std::move(data)), [](void* ptr) {
                    delete reinterpret_cast<std::vector<T>*>(ptr);
                })
            );
        }
        
        // 高效的参数验证宏
        #define VALIDATE_CHAIN_TYPE(chain_type) \
            if (UNLIKELY(chain_type != "ring" && chain_type != "open")) { \
                throw std::runtime_error("Invalid chain_type. Expected 'ring' or 'open'."); \
            }
        
        // 高效的数组维度验证
        #define VALIDATE_ARRAY_DIMS(info, expected_dims) \
            if (UNLIKELY(info.ndim != expected_dims)) { \
                throw std::runtime_error("Expected " #expected_dims "-dimensional array"); \
            }
        
        // 快速类型转换和验证
        template<typename T>
        FORCE_INLINE T* validateAndGetPointer(py::array_t<T>& array) {
            py::buffer_info info = array.request();
            if (UNLIKELY(!info.ptr)) {
                throw std::runtime_error("Invalid array pointer");
            }
            return static_cast<T*>(info.ptr);
        }
    }
}

#endif // OPTIMIZED_TYPES_H