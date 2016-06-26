EXECS=par2_diff
MPICC?=mpicc

all: ${EXECS}

par2_diff: par2_diff.c
	${MPICC} -o par2_diff par2_diff.c

clean:
	rm ${EXECS}

#mpirun -np 4 ./par2_diff

