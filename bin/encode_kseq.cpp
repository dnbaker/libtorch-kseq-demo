#include <torch/torch.h>
#include <iostream>
#include <zlib.h>
#include "kseq.h"


KSEQ_INIT(gzFile, gzread)

template <typename T>
void pretty_print(const std::string& info, T&& data) {
  std::cout << info << std::endl;
  std::cout << data << std::endl << std::endl;
}

enum Alph {
    DNA5, // N for unknown
    PROTEIN17, // 16, + end-of-sequence token
    PROTEIN22, // 21
};

struct Tokenizer {
    std::array<uint8_t, 128> lut_;
    const Alph alphabet_;
    const long long signed alphabet_size_;
    const int include_eos_, include_bos_;
    size_t full_alphabet_size() const {return alphabet_size_ + include_eos_ + include_bos_;}
    // Use 0 as eos/terminus by default
    // Even without eos, that value will often be used to pad for batches.
    // 0 = padding
    // If include eos, then padding
    Tokenizer(Alph a=DNA5, bool include_eos=true, bool include_bos=false): alphabet_(a), alphabet_size_(a == DNA5 ? 5: a == PROTEIN17 ? 16: 21), include_eos_(include_eos), include_bos_(include_bos) {
        switch(a) {
            case PROTEIN22: {
                const int baseval = 1 + include_bos_ + include_eos_;
                std::fill(lut_.begin(), lut_.end(), 0);
                static constexpr int ids [] {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 11, 4, 20, 20};
                static constexpr int chars [] {'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', 'X', 'O', 'U', 'B', 'Z'};
                for(size_t i = 0; i < sizeof(ids) / sizeof(int); ++i) {
                    lut_[chars[i]] = ids[i] + baseval;
                } 
                break;
            }
            case DNA5: {
                const int baseval = 1 + include_bos_ + include_eos_;
                std::fill(lut_.begin(), lut_.end(), 4);
                lut_['a'] = lut_['A'] = 0 + baseval;
                lut_['c'] = lut_['C'] = 1 + baseval;
                lut_['g'] = lut_['G'] = 2 + baseval;
                lut_['t'] = lut_['T'] = 3 + baseval;
                break;
            }
            default: throw std::runtime_error("Not supported");
        }
    }
    char translate(int x) const {return lut_[x];}
    int eos() const {return include_eos_;}
    int bos() const {return eos() + 1;}
    torch::Tensor tokenize(const std::string &s) {return tokenize(s.data(), s.size());}
    torch::Tensor tokenize(const char *s, int64_t m) {
        auto options = torch::TensorOptions().dtype(torch::kFloat32);
        // bos is its own symbol, but eos is the same as padding symbol
        // This way, a sequence hits EOS and then stays in EOS forever.
        torch::Tensor ret = torch::zeros({m + include_eos_ + include_bos_, alphabet_size_ + include_bos_ + include_eos_ + 1}, options);
        int64_t startind = 0;
        if(include_bos_) {
            ret.index_put_({0, int64_t(bos())}, 1);
            ++startind;
        }
        for(int64_t i = startind;i < m + startind; ++i) {
            ret.index_put_({i, int64_t(lut_[s[i - startind]])}, 1);
        }
        if(include_eos_) ret.index_put_({m + 1, eos()}, 1);
        return ret;
    }
};

int main(int argc, char **argv) {
  int alphabet = argc > 1 ? std::atoi(argv[1]): 0;
  Tokenizer tok((Alph)alphabet);
  std::string seq("ACGTACGTACGT");
  torch::Tensor seq_tokenized = tok.tokenize(seq);
  std::cerr << "seq: " << seq << " >>> " << seq_tokenized << '\n';
  // Create an eye tensor
  if(torch::cuda::is_available()) {
     std::cerr << "CUDA is available!\n" << '\n';
  } else {
     std::cerr << "No CUDA, sorry\n";
  }
  if(argc > 2) {
    int64_t maxl = 0;
    std::vector<std::pair<char *, int64_t>> seqs;
    gzFile fp = gzopen(argv[2], "rb");
    if(!fp) return 1;
    kseq_t *ks = kseq_init(fp);
    for(;kseq_read(ks) >= 0;) {
        std::pair<char *, int64_t> p((char *)std::malloc(ks->seq.l + 1), ks->seq.l);
        seqs.push_back(p);
        std::memcpy(p.first, ks->seq.s, ks->seq.l + 1);
        maxl = std::max(maxl, int64_t(ks->seq.l + 1));
    }
    gzclose(fp);
    kseq_destroy(ks);
    torch::Tensor seq_stacked = torch::zeros({int64_t(seqs.size()), maxl, int64_t(tok.alphabet_size_)});
    std::cerr << "Shape: " << seqs.size() << ", " << tok.alphabet_size_ << ", " << maxl << '\n';
    auto accessor = seq_stacked.accessor<float, 3>();
    for(size_t id = 0; id < seqs.size(); ++id) {
        for(int64_t idx = 0; idx < seqs[id].second; ++idx) {
            accessor[id][tok.translate(seqs[id].first[idx])][idx] = 1;
        }
        for(int64_t idx = seqs[id].second; idx < maxl; ++idx) {
            accessor[id][tok.eos()][idx] = 1;
        }
    }
    torch::save(seq_stacked, "stacked_seq.pt");
    for(auto &s: seqs) std::free(s.first);
  }
}
