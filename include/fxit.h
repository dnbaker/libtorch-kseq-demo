#pragma once
#include <cassert>
#include "fpwrap/fpwrap.h"
#include <iostream>
#include <functional>
#include <istream>
#include <fstream>


#ifndef unlikely
#  if __GNUC__ || __clang__
#    define unlikely(x) __builtin_expect((x), 0)
#  else
#    define unlikely(x) (x)
#  endif
#endif


namespace fx {

struct StrPair {
    std::string s_; std::string q_;
    std::string to_string() const {return s_ + '|' + q_;}
};
struct FxFile {
    std::ifstream ifs_;
    StrPair sq_;
    bool is_last_;
    std::string line_;
    char *buf_;

    struct fx_iterator;
    auto end() {return fx_iterator(*this);}
    struct fx_iterator {
        FxFile &ref_;
        fx_iterator(FxFile &ref): ref_(ref) {}
        template<typename T>
        bool operator!=(const T &x) const {
            //std::fprintf(stderr, "Containing %s. islast: %s. isgood; %s. ret: %d\n", ref_.sq_.to_string().data(), ref_.is_last_ ? "true": "false", ref_.ifs_.good() ?"true":"false", ref_.ifs_.good() || ref_.is_last_);
            return ref_.ifs_.good() || ref_.is_last_ == 0;
        }
        void increment() {
            if(ref_.line_.empty()) {
                //std::fprintf(stderr, "Ref line is empty\n");
                if(!std::getline(ref_.ifs_, ref_.line_)) {
                    assert(ref_.ifs_.eof());
                    //assert(!(*this != ref_.end()));
                    return;
                }
            }
            //std::fprintf(stderr, "Calling increment, now have line %s\n", ref_.line_.data());
            const auto c = ref_.line_[0];
            if(c == '>') {
                if(!std::getline(ref_.ifs_, ref_.seq())) throw "a party";
                while(std::getline(ref_.ifs_, ref_.line_) && ref_.line_[0] != '>') {
                    ref_.seq() += ref_.line_;
                }
                assert(ref_.line_[0] == '>' || ref_.line_[0] == '@' || ref_.ifs_.eof());
            } else if(c == '@') {
                std::getline(ref_.ifs_, ref_.seq());
                std::getline(ref_.ifs_, ref_.line_); // Skip +
                std::getline(ref_.ifs_, ref_.qual());
                std::getline(ref_.ifs_, ref_.line_);
                if(unlikely(ref_.seq().size() != ref_.qual().size())) throw std::runtime_error("Seq/Qual of differing lengths.");
                std::getline(ref_.ifs_, ref_.line_); do {ref_.seq() += ref_.line_;std::getline(ref_.ifs_, ref_.line_);} while(ref_.line_[0] != '>');
            } /*else if(c == 0) {
                assert(ref_.ifs_.eof());
            } */else {
                    throw std::runtime_error("Unexpected character!");
            }
            //std::fprintf(stderr, "Got sequence %s\n", ref_.sq_.to_string().data());
        }
        auto operator++() {
            if(ref_.ifs_.eof()) {
                //std::fprintf(stderr, "Reversing is_last from %d\n", ref_.is_last_);
                ref_.is_last_ = !ref_.is_last_;
            }
            else increment();
            return *this;
        }
        auto &operator*() {
            //if(ref_.seq().empty())
            //    increment();
            //if(ref_.ifs_.good())
            return ref_.sq_;
        }
        auto operator->(){/*if(ref_.seq().empty()) increment(); */return &ref_.sq_;}
        fx_iterator operator++(int) = delete;
    };
    std::string &seq() {return sq_.s_;}
    std::string &qual() {return sq_.q_;}
    FxFile(std::string path): ifs_(path), is_last_(false) {
        std::ios_base::sync_with_stdio(false);
        buf_ = new char[1 << 16];
        ifs_.rdbuf()->pubsetbuf(buf_, 1<<16);
    }
    ~FxFile() {delete buf_;}
    fx_iterator begin() {
        auto ret = fx_iterator(*this);
        ifs_.clear(); ifs_.seekg(0);
        return ++ret;
    }
};

} // namespace fx (just like Atlanta)
