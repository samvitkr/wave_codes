#pragma once
// git.h
// https://raw.githubusercontent.com/andrew-hardin/cmake-git-version-tracking/master/git.h
//
// Released under the MIT License.
// https://raw.githubusercontent.com/andrew-hardin/cmake-git-version-tracking/master/LICENSE

#include <stdbool.h>

#ifdef __cplusplus
#define GIT_VERSION_TRACKING_EXTERN_C_BEGIN extern "C" {
#define GIT_VERSION_TRACKING_EXTERN_C_END }
#else
#define GIT_VERSION_TRACKING_EXTERN_C_BEGIN
#define GIT_VERSION_TRACKING_EXTERN_C_END
#endif

// Don't mangle the C function names if included in a CXX file.
GIT_VERSION_TRACKING_EXTERN_C_BEGIN

/// Is the metadata populated? 
//
/// We may not have metadata if there wasn't a .git directory
/// (e.g. downloaded source code without revision history).
bool git_IsPopulated();

/// Were there any uncommitted changes that won't be reflected
/// in the CommitID?
bool git_AnyUncommittedChanges();

/// The commit SHA1.
const char* git_CommitSHA1();

/// The ISO8601 commit date.
const char* git_CommitDate();

/// The commit describe.
const char* git_Describe();

/// The symbolic reference tied to HEAD.
const char* git_Branch();

GIT_VERSION_TRACKING_EXTERN_C_END
#undef GIT_VERSION_TRACKING_EXTERN_C_BEGIN
#undef GIT_VERSION_TRACKING_EXTERN_C_END
