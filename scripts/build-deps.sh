#!/usr/bin/env bash
# Builds GMP, CLN and GiNaC from source as static libraries with position-independent code and
# installs them into PREFIX, together with their license files in PREFIX/share/licenses.
# cibuildwheel runs this in before-all (see pyproject.toml); the extension module links them statically,
# so the wheels need no system GiNaC.
#
# usage: scripts/build-deps.sh PREFIX
set -euo pipefail

GMP_VERSION=6.3.0
GMP_SHA256=a3c2b80201b89e68616f4ad30bc66aee4927c3ce50e33929ca819d5c43538898
CLN_VERSION=1.3.7
CLN_SHA256=7c7ed8474958337e4df5bb57ea5176ad0365004cbb98b621765bc4606a10d86b
GINAC_VERSION=1.8.10
GINAC_SHA256=6cac1973a5325de0b9bcb8e392988ae95fbc37aa66c0f1f1d3b8e64c08cec1b9

prefix="${1:?usage: $0 PREFIX}"
mkdir -p "${prefix}"
prefix="$(cd "${prefix}" && pwd)"
jobs="$(getconf _NPROCESSORS_ONLN 2>/dev/null || echo 2)"
work="$(mktemp -d)"
trap 'rm -rf "${work}"' EXIT

static_pic=(--prefix="${prefix}" --enable-static --disable-shared --with-pic)
export PKG_CONFIG_PATH="${prefix}/lib/pkgconfig${PKG_CONFIG_PATH:+:${PKG_CONFIG_PATH}}"

fetch() {  # fetch URL SHA256: downloads, checks and unpacks into the work directory
    local file="${work}/${1##*/}"
    curl --fail --silent --show-error --location --retry 5 --output "${file}" "$1"
    local actual
    actual="$( (sha256sum 2>/dev/null || shasum -a 256) < "${file}" | cut -d' ' -f1)"
    if [ "${actual}" != "$2" ]; then
        echo "checksum mismatch for $1: ${actual}, expected $2" >&2
        exit 1
    fi
    tar -xf "${file}" -C "${work}"
}

licenses() {  # licenses NAME DIR FILES...: copies the license files of a library
    local name="$1" dir="$2"
    shift 2
    mkdir -p "${prefix}/share/licenses/${name}"
    for f in "$@"; do cp "${dir}/${f}" "${prefix}/share/licenses/${name}/"; done
}

# GMP picks code for the CPU it is built on; on x86_64 a fat build selects it at run time instead.
gmp_flags=()
case "$(uname -m)" in
    x86_64 | amd64) gmp_flags+=(--enable-fat) ;;
esac

# GNU mirror instead of gmplib.org, which often refuses downloads from CI runners
fetch "https://ftp.gnu.org/gnu/gmp/gmp-${GMP_VERSION}.tar.xz" "${GMP_SHA256}"
(cd "${work}/gmp-${GMP_VERSION}" && ./configure "${static_pic[@]}" ${gmp_flags[@]+"${gmp_flags[@]}"} && make -j"${jobs}" && make install)
licenses gmp "${work}/gmp-${GMP_VERSION}" COPYING COPYING.LESSERv3 COPYINGv2 COPYINGv3

fetch "https://www.ginac.de/CLN/cln-${CLN_VERSION}.tar.bz2" "${CLN_SHA256}"
(cd "${work}/cln-${CLN_VERSION}" && ./configure "${static_pic[@]}" --with-gmp="${prefix}" && make -j"${jobs}" && make install)
licenses cln "${work}/cln-${CLN_VERSION}" COPYING

# compile_ex needs dlopen and a compiler at run time; the bindings do not use it.
fetch "https://www.ginac.de/ginac-${GINAC_VERSION}.tar.bz2" "${GINAC_SHA256}"
(cd "${work}/ginac-${GINAC_VERSION}" && ./configure "${static_pic[@]}" --disable-excompiler \
    && make -C ginac -j"${jobs}" && make -C ginac install && make install-pkgconfigDATA)
licenses ginac "${work}/ginac-${GINAC_VERSION}" COPYING
