#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
project_dir="$(cd "${script_dir}/.." && pwd)"
project_name="$(basename "${project_dir}")"
parent_dir="$(dirname "${project_dir}")"
default_output="${parent_dir}/${project_name}-source.zip"
output_path="${1:-${default_output}}"

case "${output_path}" in
  /*) ;;
  *) output_path="$(pwd)/${output_path}" ;;
esac

case "${output_path}" in
  "${project_dir}"/*)
    echo "error: output ZIP must be outside the project tree" >&2
    exit 2
    ;;
esac

rm -f "${output_path}"

(
  cd "${parent_dir}"
  COPYFILE_DISABLE=1 zip -q -r -X \
    "${output_path}" \
    "${project_name}" \
    -x \
    "${project_name}/build/*" \
    "${project_name}/build-*/*" \
    "${project_name}/docs/doxygen/*" \
    "${project_name}/docs/doxygen-warnings.log" \
    "${project_name}/output/*" \
    "${project_name}/.git/*" \
    "*/.DS_Store" \
    "*/._*" \
    "*/__MACOSX/*"
)

echo "Created ${output_path}"
