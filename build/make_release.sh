#!/bin/bash
#
# See 
# https://circleci.com/docs/nightly-builds/
# 
set -e
set -u
set -o pipefail

if [[ ! $CIRCLE_TOKEN ]]; then
  echo 'Set the CIRCLE_TOKEN variable in your environment'
  echo 'Check the following link to get your token https://circleci.com/account/api'
  exit 1
fi 

echo "This script trigger a T-Coffee build and makes a new release"

# Ask for confirmation 
read -u 1 -p "* Confirm which type of release you want to creare [beta/stable] " -r
if [[ $REPLY == 'beta' ]]; then
  RELEASE=0
elif [[ $REPLY == 'stable' ]]; then
  RELEASE=1
else 
  echo Unknown option. You must type 'beta' or 'stable'
  exit 1
fi 

RELEASE=1

function trigger_build() {

	local PROJECT=cbcrg/tcoffee
	local BRANCH=${2:-master}
	echo "Triggering T-Coffee build [${BRANCH}] -- RELEASE=${RELEASE}"

	local trigger_build_url=https://circleci.com/api/v1/project/${PROJECT}/tree/${BRANCH}?circle-token=${CIRCLE_TOKEN}

	local params=''
	[[ $RELEASE == 1 ]] && params+="\"RELEASE\": \"1\","
	params+="\"RUN_NIGHTLY_BUILD\": \"true\","
	params+="\"FUNCTIONAL_TEST_TARGET\": \"staging-dawn-435.herokuapp.com\""
	
	local post_data=$(cat <<EOF
	{
	  "build_parameters": { $params }
	}
	EOF)

	curl \
	--header "Accept: application/json" \
	--header "Content-Type: application/json" \
	--data "${post_data}" \
	--request POST ${trigger_build_url}
	
	echo ''
}

trigger_build