#!/usr/bin/env bash
#
# Add user to all BAM Read ACLs. This script is used to ensure that all BAMs can be read by the hmf-crunch user. It finds all the BAM files
# on which the user is not on the ACL, then adds them via the SBP API.
#
# NOTE: the environment you run this in must have the correct auth and project set throughout the run

function log() {
  echo "$(date) $@"
}

function filefilter() {
  egrep -e 'gs://hmf-.*/.*/.*/mapping/.*\.realigned\.bam(\.bai)?' -e 'gs://hmf-.*/.*/.*/aligner/.*.bam(\.bai)?' \
    -e 'gs://hmf-.*/.*/.*/mapping/.*\.realigned.*\.flagstat$' -e 'gs://hmf-.*/.*/.*/aligner/.*\.flagstat$'
}

API="https://api.hartwigmedicalfoundation.nl/hmf/v1"
CURL="curl --retry 3 -s --cert $2 --key $3"
INI_FILTER_ID=6
WORKING_DIR="$(mktemp -d)"
SHARED_FILE_IDS="${WORKING_DIR}/shared_file_ids"
SHARED_FILE_PATHS="${WORKING_DIR}/shared_file_paths"
STORED_FILES="${WORKING_DIR}/stored_files"
RUNS_FETCHED="${WORKING_DIR}/runs_fetched"
NEW_FILES="${WORKING_DIR}/new_files_to_share"
FILE_IDS_PATHS="${WORKING_DIR}/file_ids_paths"
FILE_IDS_TO_SHARE="${WORKING_DIR}/file_ids_to_share"
ALL_RUNS="${WORKING_DIR}/all_runs"

for dep in curl jq; do
  which $dep >/dev/null || (echo "Dependency $dep not found" && exit 1)
done

[[ $# -ne 4 ]] && echo "USAGE: $(basename $0) [group_id] [api.crt] [api.key] [requester-pays project name]" && exit 1
log "Fetching list of currently-shared file ids for group \"${1}\""
curl -s --cert "$2" --key "$3" "${API}/groups/${1}/files" | jq -c '.[].file_id' > $SHARED_FILE_IDS
log "Resolving file ids to paths"
for fileid in $(cat $SHARED_FILE_IDS); do
  curl -s --cert "$2" --key "$3" "${API}/files/${fileid}" | jq -r '.filepath' >> $SHARED_FILE_PATHS
done

if [[ "$(wc -l $SHARED_FILE_IDS | awk '{print $1}')" -ne "$(wc -l $SHARED_FILE_PATHS | awk '{print $1}')" ]]; then
  log "Counts do not match. Continuing but results may be incomplete."
fi

log "Fetching list of stored files"
for y in $(seq 2017 $(date '+%Y')); do
  for m in $(seq 1 12); do
    datestr="$(printf "%d-%02d" $y $m)"
    gsutil -u ${4} ls -r gs://hmf-output-${datestr}/** | filefilter >> $STORED_FILES
  done
done

log "Sorting data files"
for unsorted in $SHARED_FILE_PATHS $STORED_FILES; do
  sort --parallel=$(nproc) -u $unsorted > $unsorted.sorted
  mv ${unsorted}.sorted $unsorted
done

diff $SHARED_FILE_PATHS $STORED_FILES | awk '$0 ~ "^> .+" {print $2}' > $NEW_FILES

for status in Success Validated; do
  log "Fetching ${status} runs"
  ${CURL} "${API}/runs?status=${status}&ini_id=$INI_FILTER_ID" | jq -cr '.[].id' >> $ALL_RUNS
done

log "Fetching files for runs"
for run_id in $(cat $ALL_RUNS); do
  ${CURL} "${API}/files?run_id=${run_id}" | jq -cr '.[] | (.id|tostring) + " " + .filepath' | filefilter >> $FILE_IDS_PATHS
  [[ $? -eq 0 ]] && echo $run_id >> $RUNS_FETCHED
done

log "Mapping paths to file ids"
for new_file in $(cat $NEW_FILES); do
  awk -v nf=$new_file '$2==nf {print $1; exit 0}' < $FILE_IDS_PATHS >> $FILE_IDS_TO_SHARE
done

if [[ -s $FILE_IDS_TO_SHARE ]]; then
  log "Sharing $(wc -l $FILE_IDS_TO_SHARE) files to group"
  for new_id in $(cat $FILE_IDS_TO_SHARE); do
    ${CURL} -H "Content-Type: application/json" -X POST -d "{\"group_id\": ${1}, \"file_id\": ${new_id}}" "${API}/groups/${1}/files"  
  done
else
  log "No new files to share"
fi

log "Leaving ${WORKING_DIR} for forensics"

