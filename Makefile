all: svp.cpp svputils.o vecops.o svp.o
	g++ -o runme main.cpp svputils.cpp vecops.cpp svp.cpp
	echo "File created. './runme [your matrix]' to run."
	echo "Your matrix input should be in the format '[1,2,3] [4,5,6] [7,8,9]'."
	make test

svputils.o: svputils.cpp svputils.hpp
	g++ -c svputils.cpp

vecops.o: vecops.cpp vecops.hpp
	g++ -c vecops.cpp

svp.o: svp.cpp svp.hpp vecops.hpp svputils.hpp
	g++ -c svp.cpp

test: testbuild
	./test

testbuild: test.cpp svputils.o vecops.o svp.o
	g++ -o test test.cpp svputils.cpp vecops.cpp svp.cpp
	echo "Test file created. './test' to run."

clean:
	echo "Cleaning up the files."
	rm -f $(wildcard *.o) $(wildcard *.txt) runme test