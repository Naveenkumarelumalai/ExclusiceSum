//************************************
// Filename:Exclusive_scan
// Purpose: Implemented an exclusive scan (Prefix scan) of an array of numbers using Blelloch algorithm
//          with both multiple CPU threads as well as single thread. To do this we will use an algorithmic 
//          pattern that arises often in parallel computing: balanced trees.
// Author: Naveen Kumar Elumalai
// Date: Dec 19, 2019   
//************************************

//Header files

#include<iostream>
#include <string.h>
#include <omp.h>
#include<vector>
#include<numeric>
#include <math.h>
#include<chrono>
using namespace std;

/**
 * Funtion: Up_sweep
 * Purpose: In the Up sweep function, we traverse the tree(array) computing partial sums at internal nodes of the tree.
 * This is also known as a parallel reduction, because after this phase, the root node (the last node in the array) holds the sum of all nodes in the array.
 * 
 * Parameters:
 * @param data corresponds to the input vector
 * @param output corresponds to the partial reduced vector 
 */
vector<int> up_sweep(vector<int> data) {

	int expo = 0, expo_plus_one = 0, j = 0; //expo represents 2^i, expo_plus_one represents 2^(i+1), i and j are the loop control expressions

	vector<int> output(data);  //output vector to store the up_sweep (partial sum) value
	for (int i = 0; i < log2(data.size()); i++){  
		
		expo_plus_one = pow(2, i + 1);
		expo = pow(2, i);

		#pragma omp parallel num_threads(1) default (shared) private(j)//OPEN_MP interface construct to use multiple core for performing partial reduction
		{
		//num threads represents number of CPU that we going to use
		//default (shared) specifies that variables should be shared to all the threads by default
		//private specifies individual copy of variable should be sent to each thread
		#pragma omp for
			for (int j = 0; j < data.size() - 1; j += expo_plus_one)
				output[j + expo_plus_one - 1] = output[j + expo - 1] + output[j + expo_plus_one - 1];
		}
	}
	return output;
}

/**
 * Funtion: Down_sweep
 * Purpose: In the down sweeep function, we traverse back down the tree from the root, 
 *          using the partial sums from the reduce phase to build the scan in place on the array.
 *Parameters:
 * @param output corresponds to the partial output of downsweep function
 */
void down_sweep(vector<int> &output){

	output[output.size() - 1] = 0;//Making the last element of the array zero to start the down_sweep
	int expo = 0, expo_plus_one = 0, temp = 0,j=0;//expo represents 2^i, expo_plus_one represents 2^(i+1), i and j are the loop control expressions
	
	for (int i = log2(output.size()) - 1; i >= 0; i--){

		expo_plus_one = pow(2, i + 1);
		expo = pow(2, i);
         
		#pragma omp parallel num_threads(1) default (shared) private(temp,j)//OPEN_MP interface construct to use multiple core for down_sweep_phase
		//num threads represents number of CPU that we going to use
		//default (shared) specifies that variables should be shared to all the threads by default
		//private specifies individual copy of variable should be sent to each thread
		{			
			#pragma omp for//dynamically shares the loop control expression 'j' to individual thread
			for (j = 0; j < output.size() - 1; j += expo_plus_one){

				temp = output[j + expo - 1];
				output[j + expo - 1] = output[j + expo_plus_one - 1];
				output[j + expo_plus_one - 1] = temp + output[j + expo_plus_one - 1];
			}
		}
	}
}


/**
 * Funtion: serial_implementation
 * Purpose: Perform exclusive scan using a single CPU
 * 
 *
 * Parameters:
 * @param input corresponds to the input vector
 * @param output corresponds to the exclusive sum output
 */

void serial_implementation(vector<int>& input, vector<int>& output) {
	for (int i = 1; i < input.size(); i++)
		output[i] = input[i - 1] + output[i - 1];
}

/**
 * Funtion: Exclusive sum check function
 * Purspose: Compares the Exclusive sum result of the Parallel implementation with the serial implementation
 *
 * Parameters:
 * @param parallel_output corresponds to the parallel exclusive scan implementation using OPENMP
 * @param serial_output corresponds to the serial exclusive scan implementaion 
 * This algoritm compares the output vector of parallel implementation with the serial implementation
 */
void Exclusive_sum_check(vector<int> &parallel_output, vector<int> &serial_output){
	for (int j = 0; j < parallel_output.size(); j++){
		if (parallel_output[j] != serial_output[j]) {
			cout << "Serial and parallel implementation results doesnot match " <<endl;
			break;
		}
	}
	cout<<endl;
	cout << "Serial and parallel implementation results match " << endl;
}



//************************************
// main function
//************************************

int main() {
	
	//Total number of Elements in vector should be in 2 powers
	int size;
	cout << "Enter the number of elements in power of 2" << endl;
	cin >> size;
	vector<int> input_vector(size), serial_output_vector(size);
	
	//Filling input data in the input_vector
	for (int i = 0; i < size; i++)
		input_vector[i] = 1 + i;

	//Start the timer for timing the parallel implementation of Exclusive scan
	auto start = chrono::system_clock::now();
	
	//Parallel implementation of exclusive scan : 
	//****input vector--> up_sweep function ("Partial reduced output")--> down_sweep function(Final output vector)****
	
	//In the upsweep function, input_vector is passed and the reduced values are stored in ouput_vector
	vector<int> parallel_output_vector = up_sweep(input_vector);

	//The output_vector from up_sweep function is passed as an input to the down_sweep function to perform final reduction
	down_sweep(parallel_output_vector);
	
	//Stop the timer for parallel implementation
	auto stop = chrono::system_clock::now();

	cout << "\n";
	chrono::duration<double> elapsed_seconds = stop - start;
	cout << "Time taken in seconds for Parallel implementaion using OpenMP " << elapsed_seconds.count() << endl;

    //Start the timer for  timing the serial implementation of Exclusive scan
	start = chrono::system_clock::now();
	
	//Sequential implementation of exclusive scan 
	serial_implementation(input_vector, serial_output_vector);
	
	stop = chrono::system_clock::now();
	//Stop the timer for serial implementation
	
	cout << "\n";
	elapsed_seconds = stop - start;
	cout << "Time taken in seconds for Serial implementation " << elapsed_seconds.count() << endl;

	//Chech the reduced outputs of parallel implementaion with the serial implementaion
	Exclusive_sum_check(parallel_output_vector, serial_output_vector);
	
	return 0;
}

//************************************
// End of Program
//************************************
