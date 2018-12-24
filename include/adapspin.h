#ifndef _ADAPSIN_H__
#define _ADAPSIN_H__
#include "torch/torch.h"
#error("This doesn't work yet.")
namespace torch {

struct SpinBlock: torch::nn::Module {
    torch::nn::Linear linear;
    size_t ndim_;
    bool normalized_;
    SpinBlock(size_t input, size_t ndim=1, bool normalized=true): linear(input, input), ndim_(ndim), normalized_(normalized) {
    }
    Tensor forward(Tensor x) {
        x = linear->forward(x);
        x = x.fft(ndim_, /* normalized= */ normalized_);
        return x;
    }
    Tensor backward(Tensor x) {
        x =
        x = x.ifft(ndim_, /* normalized= */normalized_);
        return x;
    }
};

struct AdapSpin: torch::nn::Module {
    std::array<SpinBlock, 3> blocks_;
    size_t ndim;
    bool normalized_;
    AdapSpin(size_t input, size_t ndim=1, bool normalized=true):
        blocks_{SpinBlock(input, ndim, normalized), SpinBlock(input, ndim, normalized), SpinBlock(input, ndim, normalized)} {
    }
    Tensor forward(Tensor x) {
        x = blocks_[0].forward(x);
        x = blocks_[1].forward(x);
        x = blocks_[2].forward(x);
        return x;
    }
    Tensor backward(Tensor x) {
        x = blocks_[0].backward(x);
        x = blocks_[1].backward(x);
        x = blocks_[2].backward(x);
        return x;
    }
};

}

#endif /* _ADAPSIN_H__ */
