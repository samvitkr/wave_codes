#pragma once

#include <common/base/logging_fwd.h>

#include <Kokkos_Core_fwd.hpp>
#include <enum.hpp/enum.hpp>
#include <mpipp/comm.h>

#include <array>
#include <memory>
#include <variant>

// Forward declarations
namespace Kokkos {
struct LayoutLeft;
} // namespace Kokkos

namespace alps::transpose {
struct TransposeOpAssign
{
  template<typename LHS, typename RHS>
  KOKKOS_FORCEINLINE_FUNCTION static void apply(LHS& lhs, RHS const& rhs)
  {
    lhs = rhs;
  }
};

struct TransposeOpAdd
{
  template<typename LHS, typename RHS>
  KOKKOS_FORCEINLINE_FUNCTION static void apply(LHS& lhs, RHS const& rhs)
  {
    lhs += rhs;
  }
};

using TransposeOps = std::variant<TransposeOpAssign, TransposeOpAdd>;

/**
 * @brief Abstract base class for transposition.
 *
 * @tparam T The scalar type (e.g., float, double).
 * @tparam ExecSpace The Kokkos execution space type.
 *
 * @note This is an abstract class. Use create_transposer() to create
 *       concrete implementations.
 */
template<class T, class ExecSpace>
class TransposerBase
{
 protected:
  using OutType = Kokkos::View<T***, Kokkos::LayoutLeft, ExecSpace>;
  using InType  = Kokkos::View<const T***, Kokkos::LayoutLeft, ExecSpace>;

 public:
  /**
   * @brief Execute the transpose operation.
   *
   * Transposes from the input view to the output view. The number of 2D slices
   * is determined by the dimensions of the input and output views.
   *
   * @param out Output view for transposed data.
   * @param in Input view containing data to transpose.
   * @param op The transpose operation (assign or add).
   * @param space The Kokkos execution space to use.
   */
  void execute(OutType const& out,
               InType const&  in,
               TransposeOps /*tag*/,
               ExecSpace const& space) const;

  /**
   * @brief Execute the transpose operation
   *
   * Transposes `howmany` 2D slices from the input view to the output view.
   *
   * @param out Output view for transposed data.
   * @param in Input view containing data to transpose.
   * @param howmany Number of consecutive 2D slices to transpose.
   * @param op The transpose operation (assign or add).
   * @param space The Kokkos execution space to use.
   */
  void execute(OutType const&   out,
               InType const&    in,
               int              howmany,
               TransposeOps     op,
               ExecSpace const& space) const;

  virtual ~TransposerBase();

 protected:
  TransposerBase(int n0, int n1, int n_blocks);

  virtual void execute_impl(OutType const&   out,
                            InType const&    in,
                            int              howmany,
                            TransposeOps     op,
                            ExecSpace const& space) const = 0;

  int n0np, n1np;
  int np;

  Logger logger_;
};

// clang-format off
ENUM_HPP_CLASS_DECL(
  TransposeMethod,
  int,
  (Default = 0)
  (Autotune = 1)
  (Single = 10)
  (All2All = 20)
  (Point2Point = 21)
  (Point2PointSHM = 22)
)
// clang-format on

class TransposerOptions
{
 public:
  TransposeMethod method{TransposeMethod::Default};

  bool               tune_tile_sizes{false};
  std::array<int, 2> tile_sizes{0, 0};
};

template<class T, typename ExecSpace>
std::unique_ptr<TransposerBase<T, ExecSpace>>
create_transposer(const mpipp::communicator& comm,
                  int                        n0,
                  int                        n1,
                  int                        max_nz,
                  TransposerOptions          options = {});

} // namespace alps::transpose
