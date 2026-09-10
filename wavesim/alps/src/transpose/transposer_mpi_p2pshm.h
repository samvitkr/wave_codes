#pragma once

#include "transposer_base.h"

#include <mpipp/comm.h>

#include <memory>

namespace alps::transpose {

template<class T, typename ExecSpace>
class TransposerMPIPoint2PointSHM final : public TransposerBase<T, ExecSpace>
{
 private:
  using base_t = TransposerBase<T, ExecSpace>;
  using base_t::logger_;
  using base_t::n0np;
  using base_t::n1np;
  using typename base_t::InType;
  using typename base_t::OutType;

 public:
  TransposerMPIPoint2PointSHM(mpipp::communicator comm,
                              int                 n0,
                              int                 n1,
                              int                 max_nz_hint,
                              TransposerOptions   options = {});

  ~TransposerMPIPoint2PointSHM() override;

 private:
  void execute_impl(const OutType&   out,
                    const InType&    in,
                    int              howmany,
                    TransposeOps     op,
                    const ExecSpace& space) const override;

  template<typename Op>
  void execute_impl(const OutType&   out,
                    const InType&    in,
                    int              howmany,
                    const ExecSpace& space) const;

  // PIMPL - Implementation class forward declaration
  struct Impl;
  std::unique_ptr<Impl> impl;
};

} // end of namespace alps::transpose
