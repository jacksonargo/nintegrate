default: build test

test: build_test run_test

echo:
	echo hello

build:
	mkdir -p lib
	gcc -Wall -c src/nintegrate.c -Iinclude/ -o lib/nintegrate.o

build_test: build
	mkdir -p bin
	gcc -Wall -lm test/src/nintegrate_test.c lib/* -Iinclude -o bin/nintegrate_test

run_test: build_test
	bin/nintegrate_test

clean:
	rm -rf lib bin
