#!/usr/bin/env bash
set -vuo pipefail

if [[ $# -eq 1 ]] && [[ "$1" = "force" ]]; then
   rm requirements.txt dev-requirements.txt interactive-requirements.txt
fi

pip-compile --extra-index-url=http://127.0.0.1:4040 setup.cfg
pip-compile dev-requirements.in
pip-compile interactive-requirements.in
pip-sync requirements.txt dev-requirements.txt interactive-requirements.txt
