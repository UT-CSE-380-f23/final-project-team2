#!/usr/bin/env bats

executable='../src/run'

@test "verify ./run exits correctly with invalid DIM" {
#run $("${factorial_executable} -1")
run $executable heat-input-bad-dim.txt
[ "$status" -eq 1 ]
#[ "${lines[-1]}" = "You don't have a 3D solver implemented yet" ]
}

@test "verify ./run exits correctly with invalid N" {
#run $("${factorial_executable} -1")
#run $executable heat-input-bad-N.txt
run $executable ../heat-input-bad-N.txt
[ "$status" -eq 1 ]
#[ "${lines[-1]}" = "num_nodes too low!" ]
}

@test "verify ./run.sh exits correctly with invalid order" {
#run $("${factorial_executable} -1")
run $executable ../heat-input-bad-order.txt
[ "$status" -eq 1 ]
}

@test "verify ./run runs correctly in verify mode" {
#run $("${factorial_executable} -1")
run $executable heat-input-good.txt
[ "${lines[0]}" = "--> verify      = 1" ]
[ "${lines[1]}" = "--> verify      = 1" ]
[ "${lines[2]}" = "--> mode        = 0" ]
[ "${lines[3]}" = "--> N           = 100" ]
[ "${lines[4]}" = "--> DIM         = 2" ]
[ "${lines[5]}" = "--> solver      = 0" ]
[ "${lines[6]}" = "--> ORDER       = 4" ]
# skip and check the output of the last two lines
[ "${lines[-2]}" = "L2 norm of the error : 1.74651e-06" ]
[ "${lines[-1]}" = "Done outputing to file" ]
}

