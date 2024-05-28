#!/usr/bin/env bash
# -*- mode:sh; -*-

conda env export > environment.yaml

black vis.py
black main.py
