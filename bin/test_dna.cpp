#include "dataset.h"
template<typename T> class TD;
using namespace seq;

int main() {
    std::string s = "ACAGCGAGCTAGGCGGTGACCAATCAGCGGCTAATCGTCAGTCACGAAATGAGCCGCTCCCCGGCTTTATAGAGAGGCGGTACCTCCAAGATGGTAGAGTATAGGAAAGTCTCGACTTGGTGTTTAGTGCATGGCCGTGGTTATACATCC";
    at::TensorOptions options = torch::dtype(c10::kFloat);
    torch::Tensor randtense = at::randn({seq::cardinality<seq::DNA2Bit>, 150}, options);
    torch::Tensor t = at::zeros({seq::cardinality<seq::DNA2Bit>, 150}, options);
    const seq::CharTransformer ct(seq::DNA2Bit);
    for(size_t i = 0; i < s.size(); t[ct.apply(s[i])][i] = 1, ++i);
    torch::print(t);
    std::vector<std::string> svec{"ACACACTATA", "TATATATATA"};
    auto sbatch = seq::encoded_seqbatch(svec.begin(), svec.end(), /* seqlen=*/ 10);
    torch::print(sbatch);
    gzFile fp = gzopen("data/test.fna", "rb");
    if(!fp) throw 1;
    kseq_t *ks = kseq_init(fp);
    auto ktensor = seq::encoded_seqbatch(ks, 2);
    print(ktensor);
    kseq_destroy(ks);
    gzclose(fp);
    torch::nn::Linear zomg = nullptr;
}
