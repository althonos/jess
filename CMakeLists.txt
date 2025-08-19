cmake_minimum_required(VERSION 3.10)
project(Jess)

add_library(jess)
target_sources(jess PRIVATE
	src/Annulus.c
	src/Join.c
	src/Molecule.c
	src/Super.c
	src/Atom.c
	src/KdTree.c
	src/Region.c
	src/TessAtom.c
	src/Jess.c
	src/Scanner.c
	src/TessTemplate.c
	src/CandidateSet.c
	src/ResIndex.c)

add_executable(jess_bin)
target_sources(jess_bin PRIVATE src/Main.c)
target_link_libraries(jess_bin jess m)
set_property(TARGET jess_bin PROPERTY OUTPUT_NAME jess)
