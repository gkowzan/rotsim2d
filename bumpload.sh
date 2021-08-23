#!/usr/bin/env bash
bumpversion --list --allow-dirty "$1"
poetry build
poetry publish -r local
