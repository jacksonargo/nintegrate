test: 
	gcc -lm -o ./nintegrate_test nintegrate.c nintegrate_test.c
	./nintegrate_test
clean:
	rm ./nintegrate_test
