#!/usr/bin/env python3
##              _               _                    _   ____
##     ___  ___| |_ _   _ _ __ | |__   ___  __ _ ___| |_|___ \
##    / __|/ _ \ __| | | | '_ \| '_ \ / _ \/ _` / __| __| __) |
##    \__ \  __/ |_| |_| | |_) | |_) |  __/ (_| \__ \ |_ / __/
##    |___/\___|\__|\__,_| .__/|_.__/ \___|\__,_|___/\__|_____|
##                       |_|
##
## Version 0.1
##
## Usage
## -----
##   $ python3 setupbeast2.py
##   $ ./lib/beast/bin/beast --packagedir lib/packages <path/to/xml>
##
## Description
## -----------
## This script downloads and extracts the BEAST2 and the following
## packages all in the `lib` directory:
##   - remaster
##
import os
import urllib.request
import tarfile
import stat
import zipfile

# ============================================================
# Configuration
lib_dir = "lib"
packages_dir = os.path.join(lib_dir, "packages")
# List of packages to install (name, URL)
packages = {
    "remaster": "https://github.com/tgvaughan/remaster/releases/download/v2.7.2/remaster.v2.7.2.zip",
}
beast_url = "https://github.com/CompEvol/beast2/releases/download/v2.7.6/BEAST.v2.7.6.Linux.x86.tgz"
beast_archive = os.path.join(lib_dir, "BEAST.v2.7.6.Linux.x86.tgz")
# ============================================================

print("Ensuring necessary directories exist...")
if not os.path.exists(lib_dir):
    os.makedirs(lib_dir)
if not os.path.exists(packages_dir):
    os.makedirs(packages_dir)

beast_extract_dir = lib_dir
beast_bin_dir = os.path.join(lib_dir, "beast/bin")
beast_java_bin = os.path.join(lib_dir, "beast/jre/bin/java")

print("Downloading BEAST2 package...")
urllib.request.urlretrieve(beast_url, beast_archive)
print("Download complete.")

print("Extracting BEAST2 package...")
with tarfile.open(beast_archive, "r:gz") as tar:
    tar.extractall(path=beast_extract_dir)
print("Extraction complete.")

print("Setting execute permissions...")
for binary in ["beast", "beauti"]:
    os.chmod(os.path.join(beast_bin_dir, binary), stat.S_IRWXU)
os.chmod(beast_java_bin, stat.S_IRWXU)
print("Permissions set.")

for package_name, package_url in packages.items():
    print(f"Downloading package: {package_name}...")
    package_dir = os.path.join(packages_dir, package_name)
    package_archive = os.path.join(packages_dir, f"{package_name}.zip")

    if not os.path.exists(package_dir):
        os.makedirs(package_dir)

    urllib.request.urlretrieve(package_url, package_archive)
    print(f"{package_name} download complete.")

    print(f"Extracting {package_name} into dedicated directory...")
    with zipfile.ZipFile(package_archive, "r") as zip_ref:
        zip_ref.extractall(package_dir)
    print(f"{package_name} extraction complete.")

print("Setup complete.")
