#!/usr/bin/env bash
# -*- mode:sh; -*-


mkdir -p lib

wget -O lib/BEAST.v2.7.6.Linux.x86.tgz https://github.com/CompEvol/beast2/releases/download/v2.7.6/BEAST.v2.7.6.Linux.x86.tgz
tar -xzf lib/BEAST.v2.7.6.Linux.x86.tgz -C lib/

chmod 750 lib/beast/bin/beast
chmod 750 lib/beast/bin/beauti
chmod 750 lib/beast/jre/bin/java
