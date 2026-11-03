#!/usr/bin/env python3
"""Export checked values from the native CPU and CUDA analytic oracles."""
from pathlib import Path
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--repo", type=Path, required=True)
parser.add_argument("--output", type=Path, default=Path("build"))
args = parser.parse_args()
args.output.mkdir(parents=True, exist_ok=True)


def replace_once(source, old, new):
    if source.count(old) != 1:
        raise ValueError(f"expected one source anchor: {old!r}")
    return source.replace(old, new)


cpu = (args.repo / "SRC/TESTS/test_cell_generated_driver.f90").read_text()
cpu = replace_once(cpu, "only: dp => real64", "only: dp => real64, error_unit")
cpu = replace_once(cpu, "    integer :: i, j, k, order\n",
                   "    integer :: i, j, k, order\n    real(dp) :: legacy(4)\n")
cpu = replace_once(cpu, "    ! A nilpotent acceleration system",
                   "    print '(a)', 'tau,component,reference,baseline,generated,baseline_error,generated_error'\n"
                   "    ! A nilpotent acceleration system")
cpu = replace_once(cpu, "        call check('nilpotent trajectory',actual,expected)\n",
                   "        call check('nilpotent trajectory',actual,expected)\n"
                   "        legacy = z + tau*(b + matmul(amat,z))\n"
                   "        legacy = legacy + tau**2*0.5_dp*(matmul(amat,b) + matmul(matmul(amat,amat),z))\n"
                   "        do i = 1,4\n"
                   "            write (*,'(es24.16,a,i0,5(a,es24.16))') tau,',',i,',',expected(i), &\n"
                   "                ',',legacy(i),',',actual(i),',',abs(legacy(i)-expected(i)), &\n"
                   "                ',',abs(actual(i)-expected(i))\n"
                   "        end do\n")
cpu = replace_once(cpu, "    print *, 'generated cell arithmetic oracle: PASS'",
                   "    write(error_unit,*) 'generated cell arithmetic oracle: PASS'")
(args.output / "cpu_trace.f90").write_text(cpu)

cuda = (args.repo / "SRC/TESTS/test_cuda_poly2.cu").read_text()
anchor = "void close(double actual, double expected, const char *message) {\n"
cuda = replace_once(cuda, anchor, anchor +
                    '  static bool first = true;\n'
                    '  if (first) { std::puts("case,actual,reference,absolute_error"); first = false; }\n'
                    '  std::printf("\\\"%s\\\",%.17g,%.17g,%.17g\\n", message, actual, expected, std::abs(actual-expected));\n')
cuda = replace_once(cuda,
                    'std::puts("CUDA poly2 analytic, recovery, periodic and lifecycle checks passed");',
                    'std::fprintf(stderr, "CUDA poly2 analytic, recovery, periodic and lifecycle checks passed\\n");')
(args.output / "cuda_trace.cu").write_text(cuda)
