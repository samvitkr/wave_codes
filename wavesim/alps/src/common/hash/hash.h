// Copyright 2022 Peter Dimov
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#pragma once

#include <climits>
#include <cstdint>
#include <functional>
#include <utility>

namespace alps::hash_ext {

template<typename T>
struct hash;

namespace detail {
template<std::size_t Bits>
struct hash_mix_impl;

// hash_mix for 64 bit size_t
template<>
struct hash_mix_impl<64>
{
  static std::uint64_t fn(std::uint64_t x)
  {
    std::uint64_t constexpr m = (std::uint64_t(0xe9846af) << 32) + 0x9b1a615d;

    x ^= x >> 32;
    x *= m;
    x ^= x >> 32;
    x *= m;
    x ^= x >> 28;

    return x;
  }
};

// hash_mix for 32 bit size_t
template<>
struct hash_mix_impl<32>
{
  static std::uint32_t fn(std::uint32_t x)
  {
    std::uint32_t constexpr m1 = 0x21f0aaad;
    std::uint32_t constexpr m2 = 0x735a2d97;

    x ^= x >> 16;
    x *= m1;
    x ^= x >> 15;
    x *= m2;
    x ^= x >> 15;

    return x;
  }
};

inline std::size_t hash_mix(std::size_t v)
{
  return hash_mix_impl<sizeof(std::size_t) * CHAR_BIT>::fn(v);
}

} // namespace detail

template<class T>
void hash_combine(std::size_t& seed, T const& v)
{
  seed = detail::hash_mix(seed + 0x9e3779b9 + std::hash<T>()(v));
}

template<typename S, typename T>
struct hash<std::pair<S, T>>
{
  size_t operator()(const std::pair<S, T>& p) const
  {
    size_t seed = 0;
    hash_combine(seed, p.first);
    hash_combine(seed, p.second);
    return seed;
  }
};
} // namespace alps::hash_ext
