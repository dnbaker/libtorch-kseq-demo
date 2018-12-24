#include "klib/kseq.h"
#include "include/fxit.h"
#include <chrono>
#include <cstdio>
KSEQ_INIT(gzFile, gzread)

int main(int argc, char *argv[]) {
    const char * path = argc > 1 ? static_cast<const char *>(argv[1]): "hello.fx";
    size_t n = 0;
    gzFile fp = gzopen(path, "rb");
    kseq_t *ks = kseq_init(fp);
#define TP() std::chrono::high_resolution_clock::now()
    auto start = TP();
    for(int c;(c = kseq_read(ks)) >= 0;) {
        ++n;
    }
    auto end = TP();
    std::fprintf(stderr, "kseq stuff: %s, %s, %zu\n", ks->seq.s, ks->name.s, ks->seq.l);
    size_t diff = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
    std::fprintf(stderr, "Total: %zu. Time: %zu\n", n, diff);
    //fx::FxFile fx(path);
    n = 0;
    auto s2 = TP();
    for(auto &x: fx::FxFile(path)) ++n/*, std::fprintf(stderr, "strings: %s\n", x.to_string().data())*/;
    auto e2 = TP();
#if 0
    fx::FxFile::fx_iterator it = fx.begin();
    for(; it != fx.end(); ++it) {
        std::fprintf(stderr, "value: %s\n", it->to_string().data());
    }
#endif
    auto diff2 = std::chrono::duration_cast<std::chrono::nanoseconds>(e2 - s2).count();
    std::fprintf(stderr, "Total: %zu. Time: %zu\n", n, diff2);
    std::fprintf(stderr, "kseq is %lf times as fast as mine\n", static_cast<double>(diff2) / diff);
}
