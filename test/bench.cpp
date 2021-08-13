#include <benchmark/benchmark.h>
#include <cstdlib>
#include <iostream>
#include <stdint.h>
#include <unordered_map>

struct ibd_rec2_t {
    // data member
    uint32_t sid1;
    uint32_t sid2;
    uint32_t pid1;
    uint32_t pid2;
    uint8_t hid1;
    uint8_t hid2;
};

struct __attribute__((packed)) A {

    uint16_t hid2 : 1, hid1 : 1, pid2_l : 14;
    uint64_t pid2_h : 6, pid1 : 20, sid2 : 19, sid1 : 19;

    A(ibd_rec2_t &r2)
    {
        uint64_t *p = (uint64_t *) ((uint16_t *) this + 1);
        uint16_t *q = ((uint16_t *) this);

        *p = r2.sid1;
        *p <<= 19;
        *p |= r2.sid2;
        *p <<= 20;
        *p |= r2.pid1;
        *p <<= 6;
        *p |= (r2.pid2 >> 14);

        *q = (r2.pid2 & 0x3fff);
        *q <<= 2;
        *q |= ((r2.hid1 << 1) + r2.hid2);
    }
};

struct __attribute__((packed)) B {

    uint16_t hid2 : 1, hid1 : 1, pid2_l : 14;
    uint64_t pid2_h : 6, pid1 : 20, sid2 : 19, sid1 : 19;

    B(ibd_rec2_t &r2)
    {
        // this works but slow
        sid1 = r2.sid1;
        hid1 = r2.hid1;
        sid2 = r2.sid2;
        hid2 = r2.hid2;
        pid1 = r2.pid1;
        pid2_l = (r2.pid2 & 0x3fff);
        pid2_h = (r2.pid2 >> 14);
    }
};

struct __attribute__((packed)) C {

    uint16_t hid2 : 1, hid1 : 1, pid2_l : 14;
    uint64_t pid2_h : 6, pid1 : 20, sid2 : 19, sid1 : 19;

    C(ibd_rec2_t &r2)
    {

        *((uint16_t *) this) = (r2.pid2 & 0x3fff) << 2;
        *((uint16_t *) this) |= ((r2.hid1 << 1) + r2.hid2);

        *(uint64_t *) ((uint16_t *) this + 1) = r2.sid1;
        *(uint64_t *) ((uint16_t *) this + 1) <<= 19;
        *(uint64_t *) ((uint16_t *) this + 1) |= r2.sid2;
        *(uint64_t *) ((uint16_t *) this + 1) <<= 20;
        *(uint64_t *) ((uint16_t *) this + 1) |= r2.pid1;
        *(uint64_t *) ((uint16_t *) this + 1) <<= 6;
        *(uint64_t *) ((uint16_t *) this + 1) |= (r2.pid2 >> 14);
    }
};

static void
bench_pointer_on_the_fly(benchmark::State &state)
{
    // Code inside this loop is measured repeatedly
    std::vector<A> vec;
    ibd_rec2_t r2;
    for (auto _ : state) {
        vec.push_back(r2);
    }
}
// Register the function as a benchmark
BENCHMARK(bench_pointer_on_the_fly);

static void
bench_attribute_acess(benchmark::State &state)
{
    // Code inside this loop is measured repeatedly
    std::vector<B> vec;
    ibd_rec2_t r2;
    for (auto _ : state) {
        vec.push_back(r2);
    }
}
// Register the function as a benchmark
BENCHMARK(bench_attribute_acess);

static void
bench_local_pointer_variable(benchmark::State &state)
{
    // Code inside this loop is measured repeatedly
    std::vector<C> vec;
    ibd_rec2_t r2;
    for (auto _ : state) {
        vec.push_back(r2);
    }
}
// Register the function as a benchmark
BENCHMARK(bench_local_pointer_variable);

class PrepareMap : public benchmark::Fixture
{
  public:
    std::unordered_map<size_t, size_t> map;
    size_t max;
    std::vector<size_t> vec, vout;

    PrepareMap()
    {
        max = 100000;
        map.reserve(size_t(2) * max);
        for (size_t i = 0; i <= max; i++)
            map[i] = 1;
        vec.resize(10000);
        vout.reserve(10000);
        for (auto &v : vec)
            v = rand();
    }
};

class PrepareMap_non_reserve : public benchmark::Fixture
{
  public:
    std::unordered_map<size_t, size_t> map;
    size_t max;
    std::vector<size_t> vec, vout;

    void
    SetUp(const ::benchmark::State &state)
    {
        max = 100000;
        for (size_t i = 0; i <= max; i++)
            map[i] = 1;
        vec.resize(10000);
        vout.reserve(10000);
        for (auto &v : vec)
            v = rand();
    }
    void
    TearDown(const ::benchmark::State &state)
    {
    }
};

BENCHMARK_F(PrepareMap, bench_unordered_map)(benchmark::State &state)
{
    for (auto _ : state) {
        vout.clear();
        for (auto v : vec)
            vout.push_back(map[v]);
    }
}

BENCHMARK_F(PrepareMap_non_reserve, bench_unordered_map_non_reserve)
(benchmark::State &state)
{
    std::cout << "start\n";
    for (auto _ : state) {
        vout.clear();
        for (auto v : vec)
            vout.push_back(map[v]);
    }
}

BENCHMARK_MAIN();
