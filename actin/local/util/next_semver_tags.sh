#!/usr/bin/env bash

git status 
[[ $? -ne 0 ]] && echo "Not in a git checkout" && exit 1
echo

set -e

TAG_REGEX="^v?(([0-9]+)\.([0-9]+)\.([0-9]+)(-((alpha|beta)\.([0-9])+))?(_([0-9a-zA-Z-]+(\.[0-9a-zA-Z-]+)*))?)$"
latest="$(git tag -l | egrep "$TAG_REGEX" | sort -V | tail -n1)"
if [[ $latest =~ $TAG_REGEX ]]; then
    branch="$(git branch | grep '^*' | cut -d' ' -f2)"
    a=${BASH_REMATCH[2]}
    b=${BASH_REMATCH[3]}
    c=${BASH_REMATCH[4]}
    beta_version=${BASH_REMATCH[8]}
    description=${BASH_REMATCH[9]}
    final_version=${BASH_REMATCH[11]}
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
            beta_suffix="beta.1_${branch}"
        else
            beta_suffix="beta.$((beta_version+1))${description}"
        fi

        echo "Major/minor/patch: Cannot make non-beta releases on non-master branch [$branch]"
        printf "Beta branch names: "
        echo "${major}-${beta_suffix} / ${minor}-${beta_suffix} / ${patch}-${beta_suffix}"
    fi
    echo
fi
