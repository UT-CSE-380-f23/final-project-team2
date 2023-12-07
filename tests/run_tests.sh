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
[ "${lines[2]}" = "--> N           = 10" ]
[ "${lines[3]}" = "--> DIM         = 2" ]
[ "${lines[4]}" = "--> solver      = 0" ]
[ "${lines[5]}" = "--> ORDER       = 4" ]
# skip and check the output of the last two lines
#[ "${lines[-2]}" = "L2 norm of the error : 0.0499302" ]
#[ "${lines[-1]}" = "Done outputing to file" ]
}

