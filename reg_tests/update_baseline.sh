#!/bin/bash
set -x

HOMEreg="$1"
test_name="$2"
commit_num="$3"

prog_name="${HOMEreg##*/}"
base_dir="${HOMEreg}/baseline_data"
base_dir_commit="${base_dir}/${test_name}.${commit_num}"

chmod 755 "$base_dir"

if [ -d "$base_dir_commit" ];then
  chmod 777 "$base_dir_commit"
  if [ -d "$base_dir_commit/sfc" ]; then
    chmod 777 "$base_dir_commit/sfc"
  fi
  rm -fr "$base_dir_commit"
fi

mkdir -p "$base_dir_commit"

copy_and_protect() {
  for file in "$@"; do
    [ -f "$file" ] && cp "$file" "$base_dir_commit" && chmod 444 "$base_dir_commit/$file"
  done
}

case "$prog_name" in
  snow2mdl)      copy_and_protect snogrb_model ;;
  ice_blend)     copy_and_protect seaice.5min.blend ;;
  global_cycle)
    if [ "$test_name" == "c192.gsi_lndincsoilnoahmp" ]; then
      copy_and_protect *.nc gaussian_interp.*
    else
      copy_and_protect *.nc
    fi
    ;;
  grid_gen)
    copy_and_protect *.nc
    mkdir -p "$base_dir_commit/sfc"
    for file in sfc/*.nc; do
      [ -f "$file" ] && cp "$file" "$base_dir_commit/sfc" && chmod 444 "$base_dir_commit/sfc/$(basename "$file")"
    done
    chmod 555 "$base_dir_commit/sfc" ;;
  *) copy_and_protect *.nc ;;
esac

chmod 555 "$base_dir_commit"
rm -f "$base_dir/$test_name"
cd "$base_dir" && ln -fs "$test_name.$commit_num" "$test_name"
