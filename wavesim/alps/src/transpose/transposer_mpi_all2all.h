#pragma once

#include "transposer_base.h"

#include <mpipp/comm.h>

namespace alps::transpose {

template<class T, typename ExecSpace>
class TransposerMPIAll2All final : public TransposerBase<T, ExecSpace>
{
 private:
  using base_t = TransposerBase<T, ExecSpace>;
  using base_t::logger_;
  using base_t::n0np;
  using base_t::n1np;
  using typename base_t::InType;
  using typename base_t::OutType;

 public:
  TransposerMPIAll2All(mpipp::communicator comm,
                       int                 n0,
                       int                 n1,
                       int                 max_nz_hint,
                       TransposerOptions   options = {});

  ~TransposerMPIAll2All() override;

 private:
  void execute_impl(const OutType&   out,
                    const InType&    in,
                    int              howmany,
                    TransposeOps     op,
                    const ExecSpace& space) const override;

  mpipp::communicator comm_;
};

} // end of namespace alps::transpose
