#!/usr/bin/env bash
# -*- mode:sh; -*-

conda env export > environment.yaml

black visualisation.py
black main.py
