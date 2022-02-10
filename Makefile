# Makefile

MPICXX = /usr/local/bin/g++-9
CXXFLAGS = -c -Wall -std=c++11 -fopenmp

mosolid: main.o element.o gibbse.o initcomp.o massb_oxide.o solution.o CGgibbsmin.o gcell.o mo_system.o
	$(MPICXX) -Wall -fopenmp main.o element.o gibbse.o initcomp.o massb_oxide.o solution.o CGgibbsmin.o gcell.o mo_system.o -o mosolid
main.o: main.cpp
	$(MPICXX) $(CXXFLAGS) main.cpp
element.o: element.cpp
	$(MPICXX) $(CXXFLAGS) element.cpp
gibbse.o: gibbse.cpp
	$(MPICXX) $(CXXFLAGS) gibbse.cpp
initcomp.o: initcomp.cpp
	$(MPICXX) $(CXXFLAGS) initcomp.cpp
massb_oxide.o: massb_oxide.cpp
	$(MPICXX) $(CXXFLAGS) massb_oxide.cpp
solution.o: solution.cpp
	$(MPICXX) $(CXXFLAGS) solution.cpp
CGgibbsmin.o: CGgibbsmin.cpp
	$(MPICXX) $(CXXFLAGS) CGgibbsmin.cpp
gcell.o: gcell.cpp
	$(MPICXX) $(CXXFLAGS) gcell.cpp
mo_system.o: mo_system.cpp
	$(MPICXX) $(CXXFLAGS) mo_system.cpp
gcell.cpp: gcell.h
clean:
	rm -f *.o mosolid
