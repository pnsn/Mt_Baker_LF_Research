#!/bin/bash
REF_FILE='MtBaker_50km_evid_recent.csv'
echo Enter role name:
read USR

echo Enter role password:
read PW

echo Enter database name:
read DB

echo Enter search radius:
read RAD

echo Enter last EVID queried:
read LEVID

TMP_FILE='./tmp_evid_list.csv'

psql -U $USR -p $PW -d $DB -v radius="'$RAD'" -v levid="'$LEVID'" -v ocsv="'$TMP_FILE'" -f ./MtBaker_Raidus_EVID_Query.sql
qml -i $TMP_FILE -o ./MtBaker_{$RAD}km_recent.qml -S