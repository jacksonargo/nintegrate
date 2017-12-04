CC=/usr/bin/gcc
INCLUDE=include
CFILES=src/*.c
CFLAGS=-Wall -I${INCLUDE}

default: test

test: build_test run_test

echo:
	echo hello

build:
	mkdir -p lib
	${CC} ${CFLAGS} -c ${CFILES}
	mv *.o lib

build_test: build
	mkdir -p bin
	${CC} ${CFLAGS} -lm test/src/nintegrate_test.c lib/*.o -Iinclude -o bin/nintegrate_test

run_test: build_test
	bin/nintegrate_test

clean:
	rm -rf lib bin
