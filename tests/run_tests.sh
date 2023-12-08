#!/usr/bin/env bats

executable='../src/run'

# bats test_tags=tag:1
@test "verify ./run exits correctly with invalid DIM" {
#run $("${factorial_executable} -1")
run $executable ../heat-input-bad-dim.txt
[ "$status" -eq 1 ]
#[ "${lines[-1]}" = "You don't have a 3D solver implemented yet" ]
}

# bats test_tags=tag:2
@test "verify ./run exits correctly with invalid N" {
#run $("${factorial_executable} -1")
#run $executable heat-input-bad-N.txt
run $executable ../heat-input-bad-N.txt
[ "$status" -eq 1 ]
#[ "${lines[-1]}" = "num_nodes too low!" ]
}

# bats test_tags=tag:3
@test "verify ./run.sh exits correctly with invalid order" {
#run $("${factorial_executable} -1")
run $executable ../heat-input-bad-order.txt
[ "$status" -eq 1 ]
}

# bats test_tags=tag:4
@test "verify ./run runs correctly in verify mode" {
#run $("${factorial_executable} -1")
run $executable ../heat-input-good.txt
[ "${lines[0]}" = "--> verify      = 1" ]
[ "${lines[1]}" = "--> mode        = 0" ]
[ "${lines[2]}" = "--> USE_PETSC   = 0" ]
[ "${lines[3]}" = "--> N           = 10" ]
[ "${lines[4]}" = "--> DIM         = 2" ]
[ "${lines[5]}" = "--> solver      = 0" ]
[ "${lines[6]}" = "--> ORDER       = 4" ]
# skip and check the output of the last two lines
#[ "${lines[-2]}" = "L2 norm of the error : 0.0499302" ]
#[ "${lines[-1]}" = "Done outputing to file" ]
}

# this test shows that the relative difference between each element of the solution
# is less that 0.001, indicating that the results obtained using GMRES and gauss-seidel
# are similar.
# bats test_tags=tag:5
#@test "check difference between 2D, N=100, 4th order gauss-seidel and PETSc GMRES" {
#run h5diff -p 0.001 ../output/no_petsc_2d_100.h5 ../output/petsc_2d_100.h5
#[ "$status" -eq 0 ]
#}

@test "run code with and without PETSc and check difference (relative tolerance) with h5diff." {
run $executable ../heat-input-good-petsc-2-100.txt ../petsc-2-100.h5
run $executable ../heat-input-good-nopetsc-2-100.txt ../nopetsc-2-100.h5
# -p specifies to compare the relative tolerances
run h5diff -p 0.001 ../nopetsc-2-100.h5 ../petsc-2-100.h5
[ "$status" -eq 0 ]
}

@test "5hdiff with 1D, N = 10, 2nd Order gauss-seidel and GMRES." {
run $executable ../heat-input-good-petsc-1-10.txt ../petsc-1-10.h5
run $executable ../heat-input-good-nopetsc-1-10.txt ../nopetsc-1-10.h5
run h5diff -p 0.000001 ../nopetsc-1-10.h5 ../petsc-1-10.h5
[ "$status" -eq 0 ]
}

