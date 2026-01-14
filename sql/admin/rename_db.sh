#!/bin/bash

set -e

if [ $# -ne 4 ]; then
    echo "Usage: $0 <user> <host> <database> <new_database>"
    echo
    echo "Renames a MySQL database by moving all tables and recreating all views."
    echo "Set the password of the user with: export MYSQL_PWD=<password>"
    exit 1
fi

user=$1
host=$2
old_db=$3
new_db=$4

function run_sql() {
  mysql -h "$host" -u "$user" -e "$1" -sss
}

function run_sql_in_db() {
  mysql -h "$host" -u "$user" "$1" -e "$2" -sss
}

function die() {
    echo "ERROR: $1"
    exit 1
}

# Fail if new database already exists
db_exists=$(run_sql "show databases like '$new_db'")
if [ -n "$db_exists" ]; then
    die "New database already exists $new_db"
fi

# Get list of all tables
TABLES=$(run_sql "select TABLE_NAME from information_schema.tables where table_schema='$old_db' and TABLE_TYPE='BASE TABLE'")
STATUS=$?
if [ "$STATUS" != 0 ] || [ -z "$TABLES" ]; then
    die "Cannot retrieve tables from $old_db"
fi

# Create new database
echo "Creating database $new_db"
run_sql "create database $new_db"

# Move all tables from old to new database
for TABLE in $TABLES; do
    echo "Renaming table $old_db.$TABLE to $new_db.$TABLE"
   run_sql_in_db "$old_db" "SET FOREIGN_KEY_CHECKS=0; rename table $old_db.$TABLE to $new_db.$TABLE"
done

# Recreate all views in the new database (rename table does not work for views). Creation fails when a view references another view
# that hasn't moved yet. Since we cannot determine the topological order automatically, we retry to create a view in a loop until
# all have been recreated or the list does not change anymore.
views=()
while IFS= read -r view; do
  views+=("$view")
done < <(run_sql "SELECT table_name FROM information_schema.views WHERE table_schema = '$old_db'")

while [ ${#views[@]} -gt 0 ]; do
  failed=()
  for view in "${views[@]}"; do
    echo "Renaming view $old_db.$view to $new_db.$view"
    create_view=$(run_sql_in_db "$new_db" "SELECT CONCAT('CREATE OR REPLACE VIEW \`$view\` AS ', view_definition) FROM information_schema.views WHERE table_schema = '$old_db' AND table_name = '$view'" | sed 's/`'"$old_db"'`\.//g')
    if ! run_sql_in_db "$new_db" "$create_view" 2>/dev/null; then
      echo "-> failed, will retry"
      failed+=("$view")
    fi
  done

  if [ ${#failed[@]} -eq ${#views[@]} ]; then
    die "Cannot rename these views: ${views[*]}"
  fi
  views=("${failed[@]}")
done

# Drop database when indeed no more tables remain
TABLES=$(run_sql "SELECT TABLE_NAME FROM information_schema.tables WHERE table_schema='$old_db' AND TABLE_TYPE='BASE TABLE'")
if [ -z "$TABLES" ]; then
    echo "Dropping database $old_db"
    run_sql "drop database $old_db"
fi
