#pragma once
#include "torch/torch.h"
#include "kseq.h"
#include <array>
#include <cstdio>
#include <algorithm>
#if ZWRAP_USE_ZSTD
#  include "zstd_zlibwrapper.h"
#else
#  include <zlib.h>
#endif

#ifndef likely
#  if defined(__GNUC__) || defined(__clang__)
#    define likely(x) __builtin_expect((x), 1)
#  else
#    define likely(x) (x)
#  endif
#endif

#ifndef unlikely
#  if defined(__GNUC__) || defined(__clang__)
#    define unlikely(x) __builtin_expect((x), 0)
#  else
#    define unlikely(x) (x)
#  endif
#endif

namespace seq {
KSEQ_INIT(gzFile, gzread)


enum Padding: int32_t {
    LEFT_PAD  = -1,
    RIGHT_PAD =  1
};

enum Alphabet {
    DNA2Bit,
    DNA4Bit,
    Protein4Bit,
    ProteinFull // TODO: Consider stop codons.
};

// TODO:
// Add reduced AA tables
// https://github.com/seqan/seqan/tree/584ae4fbff8312cfe31b3e0aea651edfb1b580ef/include/seqan/reduced_aminoacid
// li10, Murphy10, Murphy5, Solis10, Cluster (20, 22, 24), buchfink 11
// On second thought, why not just use SeqAn?

static constexpr size_t NBITS_PER_ITEM [] {
    2,
    4,
    4,
    5
};
static constexpr size_t CARDINALITY [] {
    4,
    5,
    15,
    23
};
static constexpr uint32_t EMPTY_SYMBOL [] {
    0,
    4,
    0,
    22
};

template<Alphabet a> static constexpr size_t nbits_per_item = NBITS_PER_ITEM[a];
template<Alphabet a> static constexpr size_t items_per_byte = 8 / NBITS_PER_ITEM[a];
template<Alphabet a> static constexpr size_t cardinality;
template<>  constexpr size_t  cardinality<DNA2Bit> = 4;
template<>  constexpr size_t  cardinality<DNA4Bit> = 5;
template<>  constexpr size_t  cardinality<Protein4Bit> = 16;
template<>  constexpr size_t  cardinality<ProteinFull> = 23;

// Todo: add 3-bit encoding of AAs?

static_assert(cardinality<DNA2Bit> == 4, "Must be 4");
static_assert(cardinality<DNA4Bit> == 5, "Must be 5");
static_assert(cardinality<Protein4Bit> == 16, "Must be 16");
static_assert(cardinality<ProteinFull> == 23, "Must be 23");

// Todo: add 3-bit encoding of AAs?


using lut_t = const std::array<int8_t, 256>;

static constexpr lut_t DNA_2B_LUT [] {
#define TWO2_CONTENT \
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
TWO2_CONTENT
};
static constexpr lut_t DNA_4B_LUT [] {
#define FOURB_CONTENT \
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4
FOURB_CONTENT
};
static constexpr lut_t PROTEIN_4B_LUT [] {
#define PROT15_CONTENT \
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 3, 4, 4, 5, 6, 7, 8, 0, 9, 8, 12, 2, 0, 10, 2, 11, 2, 13, 0, 12, 14, 0, 15, 0, 0, 0, 0, 0, 0, 0, 1, 0, 3, 4, 4, 5, 6, 7, 8, 0, 9, 8, 12, 2, 0, 10, 2, 11, 2, 13, 0, 12, 14, 0, 15, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
PROT15_CONTENT
};
static constexpr lut_t PROTEIN_FULL_LUT [] {
/*
arr = [20] * 256
for ind, i in enumerate("A, B, C, D, E, F, G, H, I, K, L, M, N, P, Q, R, S, T, V, W, Y, Z".split(", ")):
    arr[ord(i)] = arr[ord(i.lower())] = ind
arr[ord('n')] = arr[ord('N')] = 20
print("    " + ", ".join(map(str, arr)))
*/
#define PROTFULL_CONTENT \
    20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 0, 1, 2, 3, 4, 5, 6, 7, 8, 20, 9, 10, 11, 22, 20, 13, 14, 15, 16, 17, 20, 18, 19, 20, 20, 21, 20, 20, 20, 20, 20, 20, 0, 1, 2, 3, 4, 5, 6, 7, 8, 20, 9, 10, 11, 22, 20, 13, 14, 15, 16, 17, 20, 18, 19, 20, 20, 21, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20
PROTFULL_CONTENT
};
static constexpr lut_t * LUTS[] {
    DNA_2B_LUT,
    DNA_4B_LUT,
    PROTEIN_4B_LUT
};
//static constexpr auto PROTEIN_4B_LUT{init_lut(Protein4Bit)};

struct StringPolicy {
    template<typename T>
    size_t get_size(const T &s) const {return s.size();}
    template<typename T>
    const char *get_string(const T &s) const {return s.data();}

    size_t get_size(const std::pair<char *, size_t> &p) const {return p.second;}
    size_t get_size(const char *s) const {return s ? std::strlen(s): size_t(0);}
    size_t get_size(const kstring_t *s) const {return s->l;}

    const char *get_string(const std::pair<char *, size_t> &p) const {return p.first;}
    const char *get_string(const char *s) const {return s;}
    const char *get_string(const kstring_t *s) const {return s->s;}
};

class CharTransformer {
protected:
    const lut_t lut_;
    const Alphabet a_;
public:
    constexpr CharTransformer(Alphabet a): lut_(a == DNA2Bit ? lut_t{TWO2_CONTENT}: a == DNA4Bit ? lut_t{FOURB_CONTENT}: lut_t{PROT15_CONTENT}), a_(a) {}
    operator const std::array<int8_t, 256> &() const {
        return lut_;
    }
    auto apply(char v) const {return lut_[static_cast<uint8_t>(v)];}
    template<typename It1, typename It2>
    void transform(It1 i1, It1 i2, It2 o1) const {
        std::transform(i1, i2, o1, [&](auto x) {return this->apply(x);});
    }
    template<typename It1>
    void transform(It1 i1, It1 i2) const {transform(i1, i2, i1);}
    /*
     Core operation for encoding.
     For generality, we use a template parameter for a StringPolicy.
     This requires providing get_size and get_string methods for the types
     dereferenced from the iterator.
    */
    template<typename It1, typename Policy=StringPolicy>
    torch::Tensor &one_hot(torch::Tensor &t, It1 i1, It1 i2, Padding pad=RIGHT_PAD, const Policy &sp=Policy()) const {
#if !NDEBUG
        for(size_t i = 0; i < t.dim(); ++i)
            std::fprintf(stderr, "size at dimension %zu is: %zu\n", i, t.size(i));
#endif
        size_t ind = 0;
        assert(i1 < i2);
        for(size_t ind = 0;i1 != i2;++ind, ++i1) {
            //print(t[ind]);
            auto &s = *i1;
            const char *str = sp.get_string(s);
            const size_t sz =  sp.get_size(s);
            std::fprintf(stderr, "Now string of the zomg: %s\n", str);
            size_t j = 0;
            if(pad == RIGHT_PAD) {
                while(j < sz)
                    t[ind][j++][apply(*str++)] = 1;
                while(j < t.size(1))      t[ind][j++][EMPTY_SYMBOL[a_]] = 1;
            } else { // LEFT_PAD
                while(j < t.size(1) - sz) t[ind][j++][EMPTY_SYMBOL[a_]] = 1;
                while(j < t.size(1))
                    t[ind][j++][apply(*str++)] = 1;
            }
        }
        return t;
    }
    template<typename It1>
    static torch::Tensor encode_onehot(Alphabet a, torch::Tensor &t, It1 i1, It1 i2, Padding pad=RIGHT_PAD) {
        CharTransformer(a).one_hot(t, i1, i2, pad);
        return t;
    }
};


template<typename It>
torch::Tensor encoded_seqbatch(It i1, It i2, int64_t seqlen=-1, Alphabet a=DNA4Bit, Padding pad=RIGHT_PAD, c10::ScalarType dt=c10::kFloat) {
    ssize_t seqlens = seqlen <= 0 ? std::max_element(i1, i2, [](const auto &x, const auto &y) {return x.size() < y.size();})->size(): seqlen;
    at::TensorOptions options = torch::dtype(dt);
    torch::Tensor ret = torch::zeros({std::distance(i1, i2), seqlens, ssize_t(CARDINALITY[a])}, options);
    print(ret);
    ret = CharTransformer::encode_onehot(a, ret, i1, i2, pad);
    return ret;
}

torch::Tensor encoded_seqbatch(kseq_t *ks, size_t batchsz, int64_t seqlen=-1, Alphabet a=DNA4Bit, Padding pad=RIGHT_PAD, c10::ScalarType dt=c10::kFloat) {
    // This reads a batch of sequences from a sequence file of $batchsz or the number of sequences in the file,
    // whichever is lesser.
    // If the sequence length is provided, then sequence of that length is used.
    // The sequence is padding on the left or right, as detailed by Padding.
    std::pair<char *, size_t> *seqs = static_cast<std::pair<char *, size_t> *>(std::malloc(batchsz * sizeof(std::pair<char *, size_t>)));
    size_t nfilled = 0;
    while(nfilled < batchsz && kseq_read(ks) >= 0) {
        seqs[nfilled++] = std::pair<char *, size_t>(ks->seq.s, ks->seq.l);
        std::fprintf(stderr, "Read seq %s\n", ks->seq.s);
        ks->seq.s = nullptr; ks->seq.l = 0; // Steal ownership, forces reallocation on next iteration.
    }
    if(seqlen < 0) seqlen = std::max_element(seqs, seqs + nfilled, [](const auto &x, const auto &y) {return x.second < y.second;})->second;
    torch::Tensor ret = torch::zeros({ssize_t(nfilled), seqlen, ssize_t(CARDINALITY[a])}, at::TensorOptions(torch::dtype(dt)));
#if !NDEBUG
    std::fprintf(stderr, "Read batch of sequences of size %zu, of requested %zu\n", nfilled, batchsz);
#endif
    ret = CharTransformer::encode_onehot(a, ret, seqs, seqs + nfilled, pad);
#if !NDEBUG
    print(ret);
#endif
    while(nfilled) std::free(seqs[--nfilled].first);
    std::free(seqs);
    return ret;
}

// TODO: consider keeping labels around
// Additionally, consider keeping track of the number of items filled each iteration.
torch::Tensor &next_seqbatch(torch::Tensor &t, std::pair<char *, size_t> *seqs, kseq_t *ks, size_t batchsz, Alphabet a=DNA4Bit, Padding pad=RIGHT_PAD, c10::ScalarType dt=c10::kFloat) {
    // This reads a batch of sequences from a sequence file of $batchsz or the number of sequences in the file,
    // whichever is lesser.
    // If the sequence length is provided, then sequence of that length is used.
    // The sequence is padding on the left or right, as detailed by Padding.
    size_t nfilled = 0;
    while(nfilled < batchsz && kseq_read(ks) >= 0) {
        if(unlikely(t.size(1) < ks->seq.l))
            throw std::runtime_error("Encountered sequence of unexpected size. Abort! (Downstream, we may allow dynamic resizing, but this requires modifications to padding and perhaps architectural changes.");
        auto &ksr = seqs[nfilled];
        if(ksr.second == ks->seq.l)
            std::memcpy(ksr.first, ks->seq.s, ksr.second);
        else if(ksr.second > ks->seq.l) std::memcpy(ksr.first, ks->seq.s, ks->seq.l), ksr.second = ks->seq.l;
        else
            ksr.first = static_cast<char *>(std::realloc(ksr.first, ks->seq.l)), std::memcpy(ksr.first, ks->seq.s, ks->seq.l), ksr.second = ks->seq.l;
        ++nfilled;
    }
#if !NDEBUG
    std::fprintf(stderr, "Read batch of sequences of size %zu, of requested %zu\n", nfilled, batchsz);
#endif
    t = CharTransformer::encode_onehot(a, t, seqs, seqs + nfilled, pad);
    return t;
}

}
