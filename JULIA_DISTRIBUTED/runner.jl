@everywhere push!(LOAD_PATH, pwd());
println("=== Compiling. ===\n");
@everywhere using main_module
# main()
main()
