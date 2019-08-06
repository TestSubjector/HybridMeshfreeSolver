@everywhere push!(LOAD_PATH, pwd());
println("=== Compiling. ===\n");
@everywhere using main_module
println("=== Compiled ===");
# main()
main()
check_leaks()