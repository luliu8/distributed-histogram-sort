#include <cmath> 
#include <algorithm>
#include <iostream>
#include <utility>
#include <queue>
#include <vector>
#include <stdio.h>
#include <cstring>
#include <assert.h>
#include <mpi.h>


#include <inttypes.h>

#include "basic_defs.h"
#include "databasics.h"
#include "solution.h"
using namespace std;

#define PROBE_MULTIPLE 24

void print_values(dist_sort_t array[], int size){
	for (int i= 0; i< size; i++){
		printf("%" PRIu64 " ",array[i]);
	}
	printf ("\n");
}


// calculate the size of data in each rank
// divide data evenly among ranks 
void get_rCounts(vector<dist_sort_size_t> &sub_counts, int nProcs,vector<dist_sort_size_t> &rCounts){
	// return what each process should have after rebalancing 
	dist_sort_size_t sum = 0;
	for (int i = 0; i < nProcs; i++){
		sum += sub_counts[i];
	}
	dist_sort_size_t base = sum / nProcs;
	for (int rank = 0; rank< nProcs; rank++){
		rCounts[rank] = rank < (sum % nProcs)? (base+1) : base;
	}
}

void rebalance(const dist_sort_t *data, const dist_sort_size_t myDataCount, dist_sort_t **rebalancedData, dist_sort_size_t *rCount) {


	// Get number of processes
	int nProcs;
	MPI_Comm_size(MPI_COMM_WORLD, &nProcs);

	// Get rank of the process
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	
	// use Allgather to get the count of elements in each processor
	vector<dist_sort_size_t> sub_counts(nProcs);
	MPI_Allgather(&myDataCount, 1, MPI_UINT64_T, sub_counts.data(), 1, MPI_UINT64_T, MPI_COMM_WORLD);
	//printf("sub_counts[rank] is %d, myDataCount is %d", sub_counts[rank], myDataCount);
	assert (sub_counts[rank] == myDataCount);
	vector<dist_sort_size_t> rCounts(nProcs);
	get_rCounts(sub_counts, nProcs, rCounts);
	*rCount = rCounts[rank]; 
	//printf ("rank %d, rCount is %d, currently has %d \n",rank, *rCount, sub_counts[rank]);
	*rebalancedData =(dist_sort_t *)malloc((*rCount) * sizeof(dist_sort_t));
	if (myDataCount == *rCount) { // dont need to send/recv anything
		//printf ("rank %d just right \n", rank);
		memcpy(*rebalancedData, data, sizeof(dist_sort_t)*myDataCount);
	} else if (myDataCount < *rCount){
		memcpy(*rebalancedData, data, sizeof(dist_sort_t)*myDataCount);
		dist_sort_size_t diff = *rCount - myDataCount;
		vector<dist_sort_t> buffer(diff);
		MPI_Status status; 
		int count; 
		while (diff > 0){
			// tag 0 is for sending data 
			MPI_Recv(buffer.data(), diff, MPI_UINT64_T, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
			MPI_Get_count(&status, MPI_UINT64_T, &count);
			//printf ("rank %d received %d \n", rank, count);
			//printf ("rank %d the last element of the received buffer %d\n",rank,  buffer[count-1]);
			//printf ("rank %d last element of rebalanced data is %d, next address is %d", rank, (*rebalancedData)[myCount-1], &((*rebalancedData)[myCount]));
			memcpy(&((*rebalancedData)[*rCount- diff]), buffer.data(), sizeof(dist_sort_t)*count);
			//printf ("rank %d , copied %d  from buffer to rebalancedData\n", rank, count);
			diff -= count; 
		}
	} else {// i have more than i need,
		memcpy(*rebalancedData, data, sizeof(dist_sort_t)*(*rCount));
		//  scan sub_counts to determine what i need to do 
		// each process use the same algorithm so they will agree on what to do 
		dist_sort_size_t send_over; 

		// process all ranks before me 
		for (int source = 0; source < rank; source++ ) {
			// if what i have is less than what i should have
			if (sub_counts[source] <= rCounts[source]) continue; 
			// if i need to send out to othr processes 
			for (int dest = 0; dest < nProcs; dest++){
				if (sub_counts[source] == rCounts[source]) break; // i'm done 
				if (sub_counts[dest] >= rCounts[dest]) continue;  // already full 
				send_over = min(sub_counts[source]- rCounts[source], rCounts[dest] - sub_counts[dest]); // decide how many to send 
				sub_counts[source] -= send_over;
				sub_counts[dest] += send_over; 
			}
		}

		//print_values(sub_counts, nProcs);
		//print_values(rCounts, nProcs);
		// figure out what i need to send 
		for (int dest = 0; dest < nProcs; dest++){
			if (sub_counts[rank] == rCounts[rank]) break; // i'm done 
			if (sub_counts[dest] >= rCounts[dest]) continue;  // already full 
			//printf ("rank %d , rCounts[rank] = %"PRIu64" sub_counts[rank]=%"PRIu64", has %"PRIu64", rank %d needs %"PRIu64" \n", rank,rCounts[rank], sub_counts[rank], rCounts[rank] - sub_counts[rank], dest, rCounts[dest] - sub_counts[dest] );
			send_over = min(sub_counts[rank]-rCounts[rank] , rCounts[dest] - sub_counts[dest]); // decide how many to send 
			MPI_Send(&data[sub_counts[rank]-send_over], (int)send_over, MPI_UINT64_T, dest, 0, MPI_COMM_WORLD);
			//printf ("rank %d sent %"PRIu64" to rank %d \n", rank , send_over, dest);

			sub_counts[rank] -= send_over; 
			//printf ("rank %d has %"PRIu64"left\n ", rank, sub_counts[rank]);
		}
		
	}

}

// counts the number of keys less or equal to each of the m probes
void get_prefix_count(dist_sort_t probes[], int num_probes, dist_sort_size_t counts[],const dist_sort_t data[]){
	dist_sort_size_t cur = 0; 
	for (int i = 0; i < num_probes; i++){
		while (data[cur] <=  probes[i]) cur++; 
		counts[i] = cur; 
	}
	return; 
}

void rank_to_binsize(dist_sort_size_t counts[], int m){
	for (int i = m-1 ; i > 0; i --){
		counts[i] -= counts[i-1];
	}
	return; 
}


// get initial guess of probes evenly spaced 
void get_initial_probes(const dist_sort_t data[], dist_sort_size_t data_size, dist_sort_t probes[], int num_probes){
	//printf ("my data : ");
	//print_values(data, 5); 
	int my_bin_size, last_probe = -1; 
	for (int i = 0; i < num_probes; i++){ 
		my_bin_size = i < data_size%(num_probes+1) ?  data_size/(num_probes+1) + 1 : data_size/(num_probes+1);
		probes[i] = data[last_probe+my_bin_size]; 
		//printf ("probe %d , %" PRIu64 " is data[%d] ", i, probes[i], last_probe+my_bin_size);
		last_probe+= my_bin_size; 
	}	
	//printf("\n");
	return;
}

/*
// get initial guess of probes evenly spaced 
void get_initial_probes(dist_sort_t probes[], int num_probes, dist_sort_t global_max){
	for (int i = 0; i < num_probes; i++){ 
		probes[i] = (i+1) * (global_max/(num_probes+1));
	}	
	return;
}
*/

void find_splitters_from_result(int num_unachieved_splitters[], dist_sort_size_t global_probe_loc[], dist_sort_t probes[], int num_probes, int num_splitters,dist_sort_t splitters[], bool splitters_achieved[], dist_sort_size_t count[], double target_bin_size){
	double t_thresh =  0.01 * target_bin_size; 
	//printf ("t_thresh %f \n", t_thresh);
	int probe_index = 0;
	for (int i = 0; i < num_splitters-1; i++){
		if (splitters_achieved[i]) continue; //already found this splitter 
		// check whether found splitter[i]
		double ideal_loc = (i+1)* target_bin_size;
		//printf ("my ideal rank %f \n", ideal_loc);
		while ((probe_index < num_probes) && (global_probe_loc[probe_index] < ideal_loc - 0.5*t_thresh)) {
			probe_index++; 
		}
		/*
		if (probe_index != num_probes){
			printf ("splitter %d ideal loc is %f,  between probes, ", i, ideal_loc);
			//printf ("%" PRIu64 " and %" PRIu64 " " ,probes[probe_index-1], probes[probe_index]);
			printf( "with global rank %" PRIu64 "and %" PRIu64 " ",global_probe_loc[probe_index-1], global_probe_loc[probe_index]);
			printf("\n");
		}*/
		if ((probe_index < num_probes) && (global_probe_loc[probe_index] < ideal_loc +  0.5*t_thresh)) {
			// found splitter
			splitters[i] = probes[probe_index];
			count[i] = global_probe_loc[probe_index]; 
			splitters_achieved[i] = true; 
			//printf ("found splitter %d ",i);
			//printf("rank %" PRIu64 " \n", i, count[i]);
		} else { 
			//case 1: probe_index == num_probes. global_probe_loc[num_probes-1] still smaller than the ideal of probes[num_probes-1]
					  // so we need to find new probes in interval num_probes==probe_index
			//case 2: probes[probe_index] is too large but probes[probe_index-1] is too small, need to find new probes in interval probe_index
			num_unachieved_splitters[probe_index]++;
		}
	}
}


// generate num probes in the interval (lower, upper]
void generate_sub_probes(dist_sort_t new_probes[] ,dist_sort_t lower, dist_sort_t upper, dist_sort_size_t num){
	dist_sort_t bin_size = (dist_sort_t)floor((upper-lower)/num);
	for (int i = 1; i< num; i++){
		new_probes[i] = lower + i * bin_size; 
	}
	new_probes[num-1] = upper; 
}



// return number of newly generated probes 
// num_unachieved_splitters: how many unachieved splitters there are in each interval 
int generate_new_probes(dist_sort_t new_probes[], int num_unachieved_splitters[], dist_sort_t old_probes[], int num_probes, int numSplitters,dist_sort_t global_max){

	int probe_index = 0; 
	int num_new_probes; 

	for (int i = 0; i< num_probes+1; i++){
		if (num_unachieved_splitters[i] == 0) continue;  // no new probes need to be generated in this interval 
		num_new_probes = PROBE_MULTIPLE * num_unachieved_splitters[i];
		if (i == 0) {
			generate_sub_probes(&new_probes[probe_index], 0, old_probes[0], num_new_probes);
		} else if (i == num_probes){
			generate_sub_probes(&new_probes[probe_index], old_probes[i-1], global_max, num_new_probes);
		} else{
			generate_sub_probes(&new_probes[probe_index], old_probes[i-1], old_probes[i], num_new_probes);
		}
		probe_index += num_new_probes; 
	}
	return probe_index; // equal to total count of new probes generated 
} 




void findSplitters(const dist_sort_t *data, const dist_sort_size_t data_size, dist_sort_t *splitters, dist_sort_size_t *counts, int numSplitters) {
	// Get number of processes
	int nProcs;
	MPI_Comm_size(MPI_COMM_WORLD, &nProcs);

	// Get rank of the process
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	//dist_sort_t my_data[data_size]; 
	//memcpy(my_data, data, sizeof(dist_sort_t)*data_size); 

	//sort(my_data, (dist_sort_size_t) data_size);

	
	//printf("rank %d , has %d  data. \n", rank, data_size);
	// can receive maximum PROBE_MULTIPLE*(numSplitters-1) probes
	int num_probes = PROBE_MULTIPLE * (numSplitters - 1);
	dist_sort_t probes[num_probes];
	dist_sort_size_t my_probe_loc[num_probes]; // the number of elements smaller or equal to each probe 
	dist_sort_t my_max = data[data_size-1];
	dist_sort_t my_min = data[0];
	//printf("rank %d  mymax is %d  mymin is %d\n", rank, my_max, my_min);  
	if (rank != 0){
		int num_received_probes;  // how many elements are received. 
		// send my data count to rank 0, so rank 0 knows global count
		MPI_Reduce(&data_size, NULL, 1, MPI_UINT64_T, MPI_SUM, 0, MPI_COMM_WORLD);
		// get my maximum and minimum value and send to rank 0 
		MPI_Reduce(&my_max, NULL, 1, MPI_UINT64_T, MPI_MAX, 0, MPI_COMM_WORLD);
		//MPI_Reduce(&my_min, NULL, 1, MPI_UINT64_T, MPI_MIN, 0, MPI_COMM_WORLD);

		// iterate until converge 
		while (true){
			MPI_Bcast(&num_received_probes, 1, MPI_INT, 0, MPI_COMM_WORLD); 
			//printf ("rank %d, received %d probes from root\n", rank, num_received_probes);
			if (num_received_probes == 0){  // converged.  no more probes to receive.
				MPI_Bcast(splitters, numSplitters, MPI_UINT64_T, 0, MPI_COMM_WORLD);
				MPI_Bcast(counts, numSplitters, MPI_UINT64_T, 0, MPI_COMM_WORLD);
				//printf ("rank %d splitters converged, first splitter is %" PRIu64 ", global bin size is %" PRIu64 "\n", rank, splitters[0], counts[0]);
				break; 
			} else { // didn't convergve, rank 0 sent probes 
				MPI_Bcast(probes, num_received_probes, MPI_UINT64_T, 0, MPI_COMM_WORLD);
				//printf ("rank %d received probes\n", rank);

				get_prefix_count(probes, num_received_probes, my_probe_loc, data);
				//printf ("rank %d my_probe_loc[0] %d, [num_received_probes-1] %d\n", rank, my_probe_loc[0], my_probe_loc[num_received_probes-1]);
				// send my probe counts to rank 0
				//printf("rank %d , my probe ranks: ", rank);
				//print_values(my_probe_loc, num_probes);
				MPI_Reduce(my_probe_loc, NULL, num_received_probes, MPI_UINT64_T, MPI_SUM, 0, MPI_COMM_WORLD);
			}
		}
	} else {
		
		// logic for rank 0 
		dist_sort_size_t global_N; 
		dist_sort_size_t global_probe_loc[num_probes];// the number of elements smaller or equal to each probe globally 
		dist_sort_t global_max;
		dist_sort_t new_probes[num_probes];
		
		bool splitters_achieved[numSplitters]= {false}; 
		
		MPI_Reduce(&data_size, &global_N, 1, MPI_UINT64_T, MPI_SUM, 0, MPI_COMM_WORLD);
		MPI_Reduce(&my_max, &global_max, 1, MPI_UINT64_T, MPI_MAX, 0, MPI_COMM_WORLD);
		//MPI_Reduce(&my_min, &global_min, 1, MPI_UINT64_T, MPI_MIN, 0, MPI_COMM_WORLD);

		splitters[numSplitters-1] = global_max;
		splitters_achieved[numSplitters-1] = true;

		//printf ("global_N %"PRIu64", global_max %"PRIu64", global_min %"PRIu64"d\n", global_N, global_max, global_min);
		double target_bin_size = 1.0 * global_N/ numSplitters; 
		//printf ("target_bin_size is %f\n", target_bin_size);

		get_initial_probes(data, data_size, probes, num_probes);
		//get_initial_probes(probes, num_probes, global_max);

		//printf ("initial probes: ");
		//print_values(probes, num_probes);
		// tell other processor how many probes i'm sending over
		MPI_Bcast(&num_probes, 1, MPI_INT, 0, MPI_COMM_WORLD); 

		while (num_probes > 0){
			MPI_Bcast(probes, num_probes, MPI_UINT64_T, 0, MPI_COMM_WORLD);
			//printf ("rank %d,broadcasted %d probes, first probe is %d\n", rank, num_probes, probes[0]);

			get_prefix_count(probes, num_probes, my_probe_loc, data);
			//printf("rank 0 , my probe ranks: ");
			//print_values(my_probe_loc, num_probes);
			//get_probe_counts_from_fellow_processes
			MPI_Reduce(my_probe_loc, global_probe_loc, num_probes, MPI_UINT64_T, MPI_SUM, 0,MPI_COMM_WORLD);
			//printf("global ranks are: ");
			//print_values(global_probe_loc, num_probes);

			/*analyize result */ 
			
			// count how many unachieved splitters are in each interval. 
			// num_unachieved_splitters[i] represents the interval (probes[i-1], probes[i]]
			// first interval is (global_min, probles[0]] , last interval is (probes[num_probes-1], global_max)
			int num_unachieved_splitters[num_probes+1] = {0}; 
			
			find_splitters_from_result(num_unachieved_splitters, global_probe_loc, probes, num_probes, numSplitters, splitters, splitters_achieved, counts, target_bin_size); 
			//printf ("finished analyzing result\n");
			// generate new probes; 
			num_probes = generate_new_probes(new_probes, num_unachieved_splitters, probes, num_probes, numSplitters, global_max); 
			//printf ("new probes: ");
			//print_values(new_probes, num_probes);
			memcpy(probes, new_probes, sizeof(dist_sort_t)*num_probes); 
			MPI_Bcast(&num_probes, 1, MPI_INT, 0, MPI_COMM_WORLD); 
		}
		// all splitters have converged at this point
		//print_values(splitters, numSplitters);


		// broadcast splitters to all processes 
		MPI_Bcast(splitters,  numSplitters, MPI_UINT64_T, 0, MPI_COMM_WORLD);
		// calculat counts and broadcast to all processes 
		// right now counts contain the rank of numSplitters-1 splitters 
		counts[numSplitters-1] = global_N; 
		rank_to_binsize(counts, numSplitters); 
		MPI_Bcast(counts, numSplitters, MPI_UINT64_T, 0, MPI_COMM_WORLD);
	}
	return; 	
}


void get_prefix_sum(int recv_counts[], int displs[], int nProcs){
	displs[0] = 0; 
	for (int i = 1; i<nProcs; i++){
		displs[i] = displs[i-1] + recv_counts[i-1];
	}
}


// assumption: numSplitters maybe upto 100x large as nProcs
//             assuming numSplitters are divisible by nProcs
void moveData(const dist_sort_t *const sendData, const dist_sort_size_t sDataCount,
		dist_sort_t **recvData, dist_sort_size_t *rDataCount,
		const dist_sort_t *const splitters, const dist_sort_t *const counts, int numSplitters) {
	// Get number of processes
	int nProcs;
	MPI_Comm_size(MPI_COMM_WORLD, &nProcs);

	// Get rank of the process
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	dist_sort_t my_data[sDataCount]; 
	memcpy(my_data, sendData, sizeof(dist_sort_t)*sDataCount);
	sort(my_data, (dist_sort_size_t) sDataCount);

	// send my data to corresponding processors
	// use 2 pointers to figure out the number of data in bin i
	int bins_per_proc = numSplitters/nProcs;

	*rDataCount = 0;
	for (int i = bins_per_proc * rank; i< bins_per_proc * (rank + 1); i++){
		*rDataCount += counts[i];
	}	

	*recvData =(dist_sort_t *)malloc((*rDataCount) * sizeof(dist_sort_t));

	int left = 0 ,right = 0; 
	/* start of all-to-all collective */
	// count how many elements need to be sent to each process
	int send_counts[nProcs];
	for (int i = 0; i < nProcs; i++){
		while (right < sDataCount && my_data[right] <= splitters[(i+1)*bins_per_proc-1]) right ++; 
		// at this point, right == sDataCount  OR my_data[right] > splitters[i] 
		// we want to send my_data[left:right] to appropriate processor
		send_counts[i] = right - left; 
		left = right; 
	}
	/*	
	if (rank == 0){
	for (int i = 0; i< nProcs; i++){
			printf("%d ", send_counts[i]);
		}
		printf("\n");
	}*/

	int recv_counts[nProcs];
	MPI_Alltoall(send_counts, 1, MPI_INT, recv_counts,1, MPI_INT, MPI_COMM_WORLD);
	
	
	int send_displ[nProcs], recv_displ[nProcs];
	get_prefix_sum(send_counts, send_displ, nProcs);
	get_prefix_sum(recv_counts, recv_displ, nProcs);


	MPI_Alltoallv(my_data, send_counts, send_displ, MPI_UINT64_T,
						*recvData, recv_counts, recv_displ,MPI_UINT64_T, MPI_COMM_WORLD);
	/* end of all-to-all collective */
	/* start of  Regular point-to-point messages 
	MPI_Request requests[nProcs]; 
	for (int i = 0; i < nProcs; i++){
		while (right < sDataCount && my_data[right] <= splitters[(i+1)*bins_per_proc-1]) right ++; 
		// at this point, right == sDataCount  OR my_data[right] > splitters[i] 
		// we want to send my_data[left:right] to appropriate processor
		int dest = i / bins_per_proc;
		MPI_Isend(&my_data[left], right - left, MPI_UINT64_T, dest, 0, MPI_COMM_WORLD, &requests[i]);
		left = right; 
	}
	int cur_index = 0, recv_count; 
	MPI_Status status;
	for (int i = 0; i < nProcs; i++){ 
		// receive data from process i 
		//printf ("rank %d, maximum can receive %"PRIu64"\n", rank, *rDataCount);
		MPI_Recv(&(*recvData)[cur_index], *rDataCount, MPI_UINT64_T, i, 0, MPI_COMM_WORLD, &status);
		MPI_Get_count(&status,MPI_UINT64_T, &recv_count);
		//printf("rank %d, received %d from rank %d\n",rank, recv_count, i);
		cur_index += recv_count; 
	}
		
	assert (cur_index == *rDataCount);

	for (int i = 0; i < nProcs; i++){
		MPI_Wait(&requests[i], MPI_STATUS_IGNORE);
	}
	end of  Regular point-to-point messages  */



	/* start of  one to all collective
	// count how many elements need to be sent to each process
	int send_counts[nProcs];
	for (int i = 0; i < nProcs; i++){
		while (right < sDataCount && my_data[right] <= splitters[(i+1)*bins_per_proc-1]) right ++; 
		// at this point, right == sDataCount  OR my_data[right] > splitters[i] 
		// we want to send my_data[left:right] to appropriate processor
		send_counts[i] = right - left; 
		left = right; 
	}

	int recv_counts[nProcs];

	for (int i = 0; i < nProcs; i++){
		// tell process i, i'm sending send_counts[i] over to you 
		MPI_Gather(&send_counts[i], 1, MPI_INT, recv_counts, 1, MPI_INT, i, MPI_COMM_WORLD);
	}

	int displs[nProcs];
	get_prefix_sum(recv_counts, displs, nProcs);

	int cur = 0; 
	for (int i = 0; i < nProcs; i++){
		MPI_Gatherv(&my_data[cur], send_counts[i], MPI_UINT64_T, *recvData, recv_counts, displs, MPI_UINT64_T, i, MPI_COMM_WORLD);
		cur += send_counts[i];
	}
	 end of  one to all collective */

	/* start of  all to all collective */
		
		


}

void sort(dist_sort_t *data, dist_sort_size_t size) {
	// You are welcome to use this sort function.
	// If you wish to implement your own, you can do that too.
	// Don't use bubblesort.

	std::sort(data,data+size);
}


