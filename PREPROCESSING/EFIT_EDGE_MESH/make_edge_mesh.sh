#!/bin/bash
# Build an EFIT mesh extended to the limiter in one step:
#   1. grid_kind 2 core (theta_geom_flux = 2) from GORILLA via dump_mesh_2d.x
#   2. edge extension with aug_edge_mesh.py
#
# usage: make_edge_mesh.sh OUTDIR N1 N3 [SFC_S_MIN] [-- aug_edge_mesh.py options]
#   N1, N3     radial / poloidal size of the core (tetra_grid.inp n1, n3)
#   SFC_S_MIN  innermost core flux surface (default 0.1). GORILLA adds extra rings
#              (n_extra_rings) between it and the first regular ring automatically.
# environment:
#   GORILLA_DIR  GORILLA source/build tree (default: the tree this script is in)
#   DUMP_X       dump_mesh_2d.x (default: $GORILLA_DIR/BUILD/dump_mesh_2d.x)
#   GFILE        g-file, relative to the run directory (default MHD_EQUILIBRIA/g_file_for_test)
#   PYTHON       python with numpy scipy matplotlib contourpy triangle (default: python3)
# example: make_edge_mesh.sh mesh_100x120 100 120 1e-3 -- --fill-axis --xref-radius 12 --plot
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
GORILLA_DIR=${GORILLA_DIR:-$(cd "$HERE/../.." && pwd)}
DUMP_X=${DUMP_X:-$GORILLA_DIR/BUILD/dump_mesh_2d.x}
GFILE=${GFILE:-MHD_EQUILIBRIA/g_file_for_test}
PYTHON=${PYTHON:-python3}

OUT=$1; N1=$2; N3=$3; shift 3
SMIN=0.1
if [ $# -gt 0 ] && [ "$1" != "--" ]; then SMIN=$1; shift; fi
[ $# -gt 0 ] && [ "$1" == "--" ] && shift

mkdir -p "$OUT/core"
OUT=$(cd "$OUT" && pwd)
RUN=$OUT/core
EX=$GORILLA_DIR/EXAMPLES/example_4
cp "$EX/field_divB0.inp" "$EX/preload_for_SYNCH.inp" "$EX/gorilla.inp" "$RUN/"
ln -sfn "$GORILLA_DIR/MHD_EQUILIBRIA" "$RUN/MHD_EQUILIBRIA"
cat > "$RUN/tetra_grid.inp" <<EOF
&TETRA_GRID_NML
 grid_kind = 2 , n1 = $N1 , n2 = 4 , n3 = $N3 ,
 sfc_s_min = $SMIN , theta_geom_flux = 2 , theta0_at_xpoint = .true. ,
 g_file_filename = '$GFILE' ,
 convex_wall_filename = 'MHD_EQUILIBRIA/convex_wall_for_test.dat' ,
 netcdf_filename = 'MHD_EQUILIBRIA/netcdf_file_for_test.nc' ,
 boole_n_field_periods = .true. , n_field_periods_manual = 1 ,
 boole_write_mesh_obj = .false. ,
 filename_mesh_rphiz = 'mesh_rphiz.obj' , filename_mesh_sthetaphi = 'mesh_sthetaphi.obj' ,
/
EOF
t0=$(date +%s.%N)
(cd "$RUN" && "$DUMP_X" > run.log 2>&1) || { echo "core dump failed, see $RUN/run.log"; exit 1; }
t1=$(date +%s.%N)
"$PYTHON" "$HERE/aug_edge_mesh.py" "$RUN" "$RUN/$GFILE" "$OUT" "$@" | tee "$OUT/generator.log"
t2=$(date +%s.%N)
printf "timing: core (GORILLA grid_kind 2 + dump) %.1f s, edge generator %.1f s\n" \
  "$(echo "$t1 - $t0" | bc)" "$(echo "$t2 - $t1" | bc)" | tee -a "$OUT/generator.log"
