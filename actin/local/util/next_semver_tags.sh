#!/usr/bin/env bash

git status 2>/dev/null
[[ $? -ne 0 ]] && echo "Not in a git checkout" && exit 1
echo

set -e

latest="$(git tag -l | egrep '^v?[0-9]+\.[0-9]+\.[0-9]+(-beta)\..+$' | sort -V | tail -n1 | sed -e 's/^v//' -e 's/-beta.*$//')"
IFS='.' read -r a b c <<< "$latest"

branch="$(git branch | grep '^*' | cut -d' ' -f2)"
if [[ $branch == "master" ]]; then 
  major="$((a+1)).0.0"
  minor="${a}.$((b+1)).0"
  patch="${a}.${b}.$((c+1))"
  beta="Cannot make beta tags on master"
else
  major="Cannot make non-beta tags on non-master branch ${branch}"
  minor="Cannot make minor tags on non-master branch ${branch}"
  patch="Cannot make patch tags on non-master branch ${branch}"
  
  local_beta="$(diff <(git tag --merged master) <(git tag --merged "$branch") | grep '^>' | sort -V | tail -n1 | cut -d' ' -f2)"
  if [[ -z $local_beta ]]; then
    beta="${a}.${b}.$((c+1))-beta.1"
  else
    local_version="${local_beta%-beta.*}"
    local_beta_version="${local_beta#${local_version}-beta\.}"
    beta="${local_version}-beta.$((local_beta_version+1))"
  fi
fi

cat<<EOM
Next tags:

  Major: $major
  Minor: $minor
  Patch: $patch
  Beta:  $beta

EOM

