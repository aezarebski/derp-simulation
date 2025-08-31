#!/usr/bin/env bash
# -*- mode:sh; -*-

black main.py

# Run the demo scripts

python src/demo-database-usage.py
python src/demo-visualise-walltimes.py
