#include<stdint.h>
#include<benchmark/benchmark.h>

struct ibd_rec2_t {
    // data member
    uint32_t sid1;
    uint32_t sid2;
    uint32_t pid1;
    uint32_t pid2;
    uint8_t hid1;
    uint8_t hid2;
};

struct __attribute__((packed)) A{

    uint16_t hid2 : 1, hid1 : 1, pid2_l : 14;
    uint64_t pid2_h : 6, pid1 : 20, sid2 : 19, sid1 : 19;

    A(ibd_rec2_t &r2)
    {
        uint64_t *p = (uint64_t *) ((uint16_t *) this + 1);
        uint16_t *q = ((uint16_t *) this);

        *p  = r2.sid1;
        *p  <<= 19;
        *p  |= r2.sid2;
        *p  <<= 20;
        *p  |= r2.pid1;
        *p  <<= 6;
        *p  |= (r2.pid2 >> 14);

        *q = (r2.pid2 & 0x3fff);
        *q <<= 2;
        *q |= ((r2.hid1 << 1) + r2.hid2);
    }
};

struct __attribute__((packed)) B{

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

struct __attribute__((packed)) C{

    uint16_t hid2 : 1, hid1 : 1, pid2_l : 14;
    uint64_t pid2_h : 6, pid1 : 20, sid2 : 19, sid1 : 19;

    C(ibd_rec2_t &r2)
    {
       
        *((uint16_t *) this) = (r2.pid2 & 0x3fff)<< 2;
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


static void bench_pointer_on_the_fly(benchmark::State& state) {
  // Code inside this loop is measured repeatedly
  std::vector<A> vec;
  ibd_rec2_t r2;
  for (auto _ : state) {
   vec.push_back(r2);
  }
}
// Register the function as a benchmark
BENCHMARK(bench_pointer_on_the_fly);

static void bench_attribute_acess(benchmark::State& state) {
  // Code inside this loop is measured repeatedly
  std::vector<B> vec;
  ibd_rec2_t r2;
  for (auto _ : state) {
   vec.push_back(r2);
  }
}
// Register the function as a benchmark
BENCHMARK(bench_attribute_acess);

static void bench_local_pointer_variable(benchmark::State& state) {
  // Code inside this loop is measured repeatedly
  std::vector<C> vec;
  ibd_rec2_t r2;
  for (auto _ : state) {
   vec.push_back(r2);
  }
}
// Register the function as a benchmark
BENCHMARK(bench_local_pointer_variable);

BENCHMARK_MAIN();
