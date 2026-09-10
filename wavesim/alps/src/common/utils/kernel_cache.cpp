#include <common/utils/kernel_cache.h>

#define XXH_INLINE_ALL
#include <xxhash.h>

#include <cassert>
#include <cstdlib>
#include <cstring>
#include <iterator>
#include <list>
#include <memory>
#include <mutex>
#include <optional>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>

namespace {

void update_xxh128(XXH3_state_t* state, void const* data, std::size_t size)
{
  assert(state != nullptr);
  XXH_errorcode const err = XXH3_128bits_update(state, data, size);
  assert(err == XXH_OK && "XXH3_128bits_update failed");
  (void)err;
}

void update_size(XXH3_state_t* state, std::size_t value)
{
  update_xxh128(state, &value, sizeof(value));
}

} // namespace

namespace alps::utils {

#ifndef NDEBUG
struct KernelCacheKey
{
  std::string              source;
  std::vector<std::string> options;

  bool operator==(KernelCacheKey const& other) const noexcept
  {
    return source == other.source && options == other.options;
  }
};

static KernelCacheKey
make_kernel_cache_key(std::string_view                     source,
                      std::vector<std::string_view> const& options)
{
  KernelCacheKey key;
  key.source = source;
  key.options.reserve(options.size());
  for (std::string_view const opt : options) {
    key.options.emplace_back(opt);
  }
  return key;
}
#endif

class KernelCacheStoreImpl
{
 public:
  using HashedKey = XXH128_hash_t;

  struct HashedKeyHash
  {
    std::size_t operator()(HashedKey const& key) const noexcept
    {
      auto const lo = key.low64;
      auto const hi = key.high64;
      return static_cast<std::size_t>(
        lo ^ (hi + 0x9e3779b97f4a7c15ULL + (lo << 6U) + (lo >> 2U)));
    }
  };

  struct HashedKeyEqual
  {
    bool operator()(HashedKey const& lhs, HashedKey const& rhs) const noexcept
    {
      return XXH128_isEqual(lhs, rhs) == 1;
    }
  };

  struct Entry
  {
    HashedKey key;
#ifndef NDEBUG
    KernelCacheKey debug_key;
#endif
    std::vector<char> binary;
  };

  using LruList = std::list<Entry>;
  using Map     = std::
    unordered_map<HashedKey, LruList::iterator, HashedKeyHash, HashedKeyEqual>;

  explicit KernelCacheStoreImpl(std::size_t capacity)
    : capacity_(capacity)
    , hash_state_(XXH3_createState())
  {
    if (hash_state_ == nullptr) {
      throw std::bad_alloc();
    }
  }

  ~KernelCacheStoreImpl() { release_hash_state(); }

  KernelCacheStoreImpl(KernelCacheStoreImpl const&)            = delete;
  KernelCacheStoreImpl& operator=(KernelCacheStoreImpl const&) = delete;

  KernelCacheStoreImpl(KernelCacheStoreImpl&& other) noexcept
    : capacity_(other.capacity_)
    , hash_state_(other.hash_state_)
    , lru_list_(std::move(other.lru_list_))
  {
    rebuild_map();
    other.map_.clear();
    other.lru_list_.clear();
    other.capacity_   = 0;
    other.hash_state_ = nullptr;
  }

  KernelCacheStoreImpl& operator=(KernelCacheStoreImpl&& other) noexcept
  {
    if (this == &other) {
      return *this;
    }

    std::scoped_lock lock(mutex_, other.mutex_);
    release_hash_state();

    capacity_   = other.capacity_;
    hash_state_ = other.hash_state_;
    lru_list_   = std::move(other.lru_list_);
    rebuild_map();

    other.map_.clear();
    other.lru_list_.clear();
    other.capacity_   = 0;
    other.hash_state_ = nullptr;
    return *this;
  }

  HashedKey make_hashed_key(std::string_view                     source,
                            std::vector<std::string_view> const& options)
  {
    assert(hash_state_ != nullptr);

    XXH_errorcode const reset_err = XXH3_128bits_reset(hash_state_);
    assert(reset_err == XXH_OK && "XXH3_128bits_reset failed");
    (void)reset_err;

    update_size(hash_state_, source.size());
    update_xxh128(hash_state_, source.data(), source.size());

    update_size(hash_state_, options.size());
    for (auto const& option : options) {
      update_size(hash_state_, option.size());
      update_xxh128(hash_state_, option.data(), option.size());
    }

    return XXH3_128bits_digest(hash_state_);
  }

  void resize(std::size_t capacity)
  {
    capacity_ = capacity;
    evict_to_capacity();
  }

  void evict_to_capacity()
  {
    while (map_.size() > capacity_) {
      assert(!lru_list_.empty()
             && "cache and map are out of sync during eviction");
      auto const lru_it = std::prev(lru_list_.end());
#ifndef NDEBUG
      std::size_t const erased = map_.erase(lru_it->key);
      assert(erased == 1
             && "LRU key not found in map_: cache structures are out of sync");
#else
      map_.erase(lru_it->key);
#endif
      lru_list_.erase(lru_it);
    }
  }

  std::size_t        capacity_{128};
  mutable std::mutex mutex_;
  XXH3_state_t*      hash_state_{nullptr};
  LruList            lru_list_;
  Map                map_;

 private:
  void rebuild_map()
  {
    map_.clear();
    for (auto it = lru_list_.begin(); it != lru_list_.end(); ++it) {
      auto const [map_it, ok] = map_.emplace(it->key, it);
      (void)map_it;
      assert(ok && "duplicate cache key encountered while rebuilding map_");
    }
    evict_to_capacity();
  }

  void release_hash_state() noexcept
  {
    if (hash_state_ != nullptr) {
      XXH_errorcode const free_err = XXH3_freeState(hash_state_);
      assert(free_err == XXH_OK && "XXH3_freeState failed");
      (void)free_err;
      hash_state_ = nullptr;
    }
  }
};

KernelCacheStore::KernelCacheStore(std::size_t capacity)
  : impl_(std::make_unique<KernelCacheStoreImpl>(capacity))
{}

KernelCacheStore::~KernelCacheStore() = default;

void KernelCacheStore::resize(std::size_t capacity)
{
  std::lock_guard<std::mutex> lock(impl_->mutex_);
  impl_->resize(capacity);
}

std::size_t KernelCacheStore::capacity() const noexcept
{
  std::lock_guard<std::mutex> lock(impl_->mutex_);
  return impl_->capacity_;
}

std::optional<KernelCacheStore::Binary>
KernelCacheStore::lookup_binary(std::string_view                     source,
                                std::vector<std::string_view> const& options)
{
  std::lock_guard<std::mutex> lock(impl_->mutex_);
  auto const hashed_key = impl_->make_hashed_key(source, options);
  auto       it         = impl_->map_.find(hashed_key);
  if (it == impl_->map_.end()) {
    return std::nullopt;
  }

#ifndef NDEBUG
  assert(it->second->debug_key == make_kernel_cache_key(source, options)
         && "KernelCacheKey hash collision detected");
#endif

  auto const& binary = it->second->binary;
  // An empty binary is treated as a miss to avoid passing an invalid module
  // image to downstream loaders.
  if (binary.empty()) {
    return std::nullopt;
  }

  // NOLINTNEXTLINE(cppcoreguidelines-no-malloc)
  auto* raw_binary_copy = static_cast<char*>(std::malloc(binary.size()));
  if (raw_binary_copy == nullptr) {
    return std::nullopt;
  }

  std::memcpy(raw_binary_copy, binary.data(), binary.size());

  assert(it->second != impl_->lru_list_.end()
         && "LRU iterator not found for cache hit");
  impl_->lru_list_.splice(
    impl_->lru_list_.begin(), impl_->lru_list_, it->second);

  return Binary{Binary::DataPtr(raw_binary_copy), binary.size()};
}

void KernelCacheStore::store(std::string_view                     source,
                             std::vector<std::string_view> const& options,
                             std::vector<char>                    binary)
{
  std::lock_guard<std::mutex> lock(impl_->mutex_);
  if (impl_->capacity_ == 0) {
    return;
  }

  auto const hashed_key = impl_->make_hashed_key(source, options);
  auto       it         = impl_->map_.find(hashed_key);
#ifndef NDEBUG
  auto key = make_kernel_cache_key(source, options);
#endif
  if (it != impl_->map_.end()) {
#ifndef NDEBUG
    assert(it->second->debug_key == key
           && "KernelCacheKey hash collision detected");
#endif
    it->second->binary = std::move(binary);
#ifndef NDEBUG
    it->second->debug_key = std::move(key);
#endif
    assert(it->second != impl_->lru_list_.end()
           && "LRU iterator not found for cache update");
    impl_->lru_list_.splice(
      impl_->lru_list_.begin(), impl_->lru_list_, it->second);
    return;
  }

  KernelCacheStoreImpl::Entry entry{hashed_key,
#ifndef NDEBUG
                                    std::move(key),
#endif
                                    std::move(binary)};
  impl_->lru_list_.emplace_front(std::move(entry));
  impl_->map_[hashed_key] = impl_->lru_list_.begin();

  impl_->evict_to_capacity();
}

void KernelCacheStore::clear()
{
  std::lock_guard<std::mutex> lock(impl_->mutex_);
  impl_->map_.clear();
  impl_->lru_list_.clear();
}

} // namespace alps::utils
