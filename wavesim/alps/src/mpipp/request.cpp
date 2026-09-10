//
// Created by xuanx004 on 3/15/24.
//

#include "request.h"

#include <mpipp/environment.h>
#include <mpipp/utility.h>

#include <spdlog/spdlog.h>

namespace mpipp {

irequest::~irequest()
{
  if (req != MPI_REQUEST_NULL) {
    logger()->warn("Destroying request without waiting or cancelling it");
    MPI_Request_free(&req);
  }
}

void irequest::cancel()
{
  MPI_Cancel(&req);
  post_hook = nullptr;
}

std::optional<status> irequest::test()
{
  int    result{};
  status s;
  MPI_Test(&req, &result, &s);

  if (result == 0) return std::nullopt;
  if (post_hook) {
    post_hook();
    post_hook = nullptr;
  }
  return s;
}

status irequest::wait()
{
  status s;
  MPI_Wait(&req, &s);
  if (post_hook) {
    post_hook();
    post_hook = nullptr;
  }
  return s;
}

std::optional<status> irequest::get_status() const
{
  int    result{};
  status s;
  MPI_Request_get_status(req, &result, &s);
  if (result == 0) return std::nullopt;
  return s;
}

irequest_pool::~irequest_pool()
{
  bool any_active{false};
  for (auto& req : reqs) {
    if (req != MPI_REQUEST_NULL) {
      any_active = true;
      MPI_Request_free(&req);
    }
  }
  if (any_active) {
    logger()->warn("Destroying request pool with active requests");
  }
}

void irequest_pool::reserve(irequest_pool::size_type count)
{
  reqs.reserve(count);
  stats.reserve(count);
  post_hooks.reserve(count);
}

void irequest_pool::cancel(size_type i)
{
  MPI_Cancel(&(reqs.at(i)));
  if (post_hooks.at(i)) {
    post_hooks[i] = nullptr;
  }
}

void irequest_pool::cancelall()
{
  for (size_type i = 0; i < reqs.size(); ++i) {
    cancel(i);
  }
}

void irequest_pool::push(irequest_pool&& other)
{
  auto old_size = size();
  auto new_size = old_size + other.size();
  reqs.resize(new_size);
  stats.resize(new_size);
  post_hooks.resize(new_size);
  for (size_type i = 0; i < other.size(); ++i) {
    reqs[old_size + i] = other.reqs[i];
    stats[old_size + i] = other.stats[i];
    post_hooks[old_size + i] = std::move(other.post_hooks[i]);
  }
  other.reqs.clear();
  other.stats.clear();
  other.post_hooks.clear();
}

void irequest_pool::push(irequest&& other)
{
  reqs.push_back(other.req);
  other.req = MPI_REQUEST_NULL;
  stats.emplace_back();
  post_hooks.push_back(std::move(other.post_hook));
  other.post_hook = nullptr;
}

void irequest_pool::push(irequest&& other, std::function<void()> post_hook)
{
  reqs.push_back(other.req);
  other.req = MPI_REQUEST_NULL;
  stats.emplace_back();
  other.post_hook = nullptr;
  post_hooks.push_back(std::move(post_hook));
}

void irequest_pool::push(MPI_Request req, std::function<void()> post_hook)
{
  reqs.push_back(req);
  stats.emplace_back();
  post_hooks.push_back(std::move(post_hook));
}

std::optional<irequest_pool::size_type> irequest_pool::waitany()
{
  int    index{};
  status s;
  MPI_Waitany(to_int_size(size()), reqs.data(), &index, &s);

  if (index == MPI_UNDEFINED) return std::nullopt;
  stats[index] = s;
  if (post_hooks[index]) {
    post_hooks[index]();
    post_hooks[index] = nullptr;
  }
  return static_cast<size_type>(index);
}

std::optional<irequest_pool::size_type> irequest_pool::testany()
{
  int    index{}, flag{};
  status s;
  MPI_Testany(to_int_size(size()), reqs.data(), &index, &flag, &s);
  if (flag != 0 and index != MPI_UNDEFINED) {
    stats[index] = s;
    if (post_hooks[index]) {
      post_hooks[index]();
      post_hooks[index] = nullptr;
    }
    return static_cast<size_type>(index);
  }
  return std::nullopt;
}

void irequest_pool::waitall()
{
  MPI_Waitall(to_int_size(size()), reqs.data(), stats.data());
  for (auto& post_hook : post_hooks) {
    if (post_hook) {
      post_hook();
      post_hook = nullptr;
    }
  }
}

bool irequest_pool::testall()
{
  int flag{};
  MPI_Testall(to_int_size(size()), reqs.data(), &flag, stats.data());
  if (flag != 0) {
    for (auto& post_hook : post_hooks) {
      if (post_hook) {
        post_hook();
        post_hook = nullptr;
      }
    }
  }
  return flag != 0;
}

std::vector<int> irequest_pool::waitsome()
{
  std::vector<int>    indices(size());
  std::vector<status> ss(size());
  int                 count;
  MPI_Waitsome(
    to_int_size(size()), reqs.data(), &count, indices.data(), ss.data());
  if (count != MPI_UNDEFINED) {
    for (int i = 0; i < count; ++i) {
      stats[indices[i]] = ss[i];
      if (post_hooks[indices[i]]) {
        post_hooks[indices[i]]();
        post_hooks[indices[i]] = nullptr;
      }
    }
    indices.resize(count);
    return indices;
  }
  return {};
}

std::vector<int> irequest_pool::testsome()
{
  std::vector<int>    indices(size());
  std::vector<status> ss(size());
  int                 count;
  MPI_Testsome(
    to_int_size(size()), reqs.data(), &count, indices.data(), ss.data());
  if (count != MPI_UNDEFINED) {
    for (int i = 0; i < count; ++i) {
      stats[indices[i]] = ss[i];
      if (post_hooks[indices[i]]) {
        post_hooks[indices[i]]();
        post_hooks[indices[i]] = nullptr;
      }
    }
    indices.resize(count);
    return indices;
  }
  return {};
}
} // namespace mpipp
