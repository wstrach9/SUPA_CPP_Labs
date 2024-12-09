CC=g++ #Name of compiler
FLAGS=-std=c++20 -w #Compiler flags (the s makes it silent)
TARGET=assignment2 #Executable name
OBJECTS=assignment2.o FiniteFunctions.o
LIBS=-I ../../GNUplot/ -lboost_iostreams

#First target in Makefile is default

${TARGET}:${OBJECTS} #Make target from objects
	@echo "Linking..."
	${CC} ${FLAGS} ${OBJECTS} ${LIBS} -o ${TARGET}
	@make clean

assignment2.o : assignment2.cxx
	@echo "Linking..."
	${CC} ${FLAGS} ${LIBS} -c assignment2.cxx


FiniteFunctions.o : FiniteFunctions.cxx
	${CC} ${FLAGS} ${LIBS} -c FiniteFunctions.cxx


FiniteFunctions.cxx :
	touch FiniteFunctions.cxx


FiniteFunctions.h :
	touch FiniteFunctions.h


assignment2.cxx :
	touch assignment2.cxx


clean: #No targets just run shell command to rm object files and emacs backups
	@rm -f *.o *~

cleantarget: #Delete the exectuable
	@rm -f ${TARGET}
