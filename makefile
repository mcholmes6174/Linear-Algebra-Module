OBJECTS = ConjugateGrad.o GaussianElim.o GMRES.o InOut.o LinAlgDriver.o \
          LinAlgToolkit.o LSQ.o Matrix.o Vector.o

.PHONY: clean

LinAlg.exe: $(OBJECTS)
	g++ $(OBJECTS) -o LinAlg.exe
# for debugging
# g++ -pedantic-errors -Wall -Weffc++ -Wextra -Wsign-conversion -std=c++2a \
#     $(OBJECTS) -o LinAlg.exe

%.o: %.cpp
	g++ -c -std=c++2a $<
# for debugging
# g++ -c -pedantic-errors -Wall -Weffc++ -Wextra -Wsign-conversion -std=c++2a $<

clean:
	rm *.o *.exe *.dat
