#!/usr/bin/env bash

git status 
[[ $? -ne 0 ]] && echo "Not in a git checkout" && exit 1
echo

set -e

latest="$(git tag -l | egrep '^v?[0-9]+\.[0-9]+\.[0-9]+(-beta\..+)?$' | sort -V | tail -n1 | sed -e 's/^v//' -e 's/-beta.*$//')"
IFS='.' read -r a b c <<< "$latest"
a=${a:-0}; b=${b:-0}; c=${c:-0}

branch="$(git branch | grep '^*' | cut -d' ' -f2)"
major="$((a+1)).0.0"
minor="${a}.$((b+1)).0"
patch="${a}.${b}.$((c+1))"

printf "Next tags:\n\n"
if [[ $branch == "master" ]]; then
  echo "Major/minor/patch: $major / $minor / $patch"
  echo "Beta branch names: Cannot make beta releases on [master]"
else
  local_beta="$(diff <(git tag --merged master) <(git tag --merged "$branch") | grep '^>' | sort -V | tail -n1 | cut -d' ' -f2)"
  if [[ -z $local_beta ]]; then
    beta_suffix="beta.1"
  else
    local_version="${local_beta%-beta.*}"
    local_beta_version="${local_beta#${local_version}-beta\.}"
    beta_suffix="beta.$((local_beta_version+1))"
    next_local_beta_version="${local_version}-${beta_suffix}"
  fi
  echo "Major/minor/patch: Cannot make non-beta releases on non-master branch [$branch]"
  printf "Beta branch names: "
  if [[ -z $next_local_beta_version ]]; then
    echo "${major}-${beta_suffix} / ${minor}-${beta_suffix} / ${patch}-${beta_suffix}"
  else
    echo "${next_local_beta_version}"
  fi
fi
echo
