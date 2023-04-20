#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h> // required for memset
#include <mpi.h>

int main(int argc, char *argv[]) {
    int size, rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	// run program with: ./ex3 MESSAGE_SIZE_BYTES
	int msg_size = 1;
	if (argc >= 2) 
		msg_size = atoi(argv[1]);

	if (rank == 0) {
		char *message = malloc(sizeof(char) * msg_size);
		memset(message, 1, msg_size);
		
		double start_time = MPI_Wtime();

		MPI_Send(message, msg_size, MPI_BYTE, 1, 1, MPI_COMM_WORLD);
		MPI_Recv(message, msg_size, MPI_BYTE, 1, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

		double stop_time = MPI_Wtime();

		printf("Time taken to transfer and recieve %d bytes is %lf\n", msg_size, (stop_time - start_time));

		// quick checksum (since array recieved is all 0, should be 0)
		int sum = 0;
		for (int i = 0; i < msg_size; i++) {
			sum += message[i];
		}
		if (sum != 0) fprintf(stderr, "There has been an error with the message recieved on rank 0.\n");
	} else if (rank == 1) {
		char *recv_buffer = malloc(sizeof(char) * msg_size);
		char *message = malloc(sizeof(char) * msg_size);
		memset(message, 0, msg_size);
		
		MPI_Recv(recv_buffer, msg_size, MPI_BYTE, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Send(message, msg_size, MPI_BYTE, 0, 2, MPI_COMM_WORLD);

		// quick checksum (since array recieved is all 1, should be msg_size)
		int sum = 0;
		for (int i = 0; i < msg_size; i++) {
			sum += recv_buffer[i];
		}
		if (sum != msg_size) fprintf(stderr, "There has been an error with the message recieved on rank 1.\n");
	} else {
		// Since we're only checking point-to-point traffic, do nothing on ranks above 1
	}

    MPI_Finalize();
}
