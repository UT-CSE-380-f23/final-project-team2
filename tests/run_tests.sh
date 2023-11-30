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
run $executable ../heat-input-good.txt
# skip and check the output of the last two lines
[ "${lines[-2]}" = "L2 norm of the error : 1.74651e-06" ]
[ "${lines[-1]}" = "Done outputing to file" ]
}

