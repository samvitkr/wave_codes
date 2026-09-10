#pragma once

#include "transposer_base.h"

#include <array>

namespace alps::transpose {

template<class T, typename ExecSpace>
class TransposerSingle final : public TransposerBase<T, ExecSpace>
{
 private:
  using base_t = TransposerBase<T, ExecSpace>;
  using base_t::logger_;
  using base_t::n0np;
  using base_t::n1np;
  using typename base_t::InType;
  using typename base_t::OutType;

 public:
  TransposerSingle(int               n0,
                   int               n1,
                   int               max_nz_hint,
                   TransposerOptions options = {});

  ~TransposerSingle() override = default;

 private:
  void execute_impl(const OutType&   out,
                    const InType&    in,
                    int              howmany,
                    TransposeOps     op,
                    const ExecSpace& space) const override;

  std::array<int, 2> tile_sizes_{16, 16};
};
} // end of namespace alps::transpose
