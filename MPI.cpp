#include <vector>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <numeric>
#include <algorithm>
#include <stdexcept>
#include <mpi.h>
#include <string>

void bruteForce(std::vector < std::vector < double > > positions, std::vector<int>& neighbors, int N, int i, double r);
void bruteForceWithoutMPI(std::vector < std::vector < double > > positions, int N, double L, double r);
double calc_square_distance(std::vector<double> pos_1, std::vector<double> pos_2);
std::vector < std::vector < double > > read_xyz_file(std::string filename, int& N, double& L);
void getMinMaxAvg(std::vector<int>& arr , int N, int &max, int &min, double &avg);
void bruteForceWithMPI(std::vector < std::vector < double > > positions, std::vector<int>& neighbors, int procId, int numProcs, int N, double r);

// int main(int argc, char **argv)
// {
//   // initialise MPI
//   MPI_Init(&argc, &argv);

//   //Choose r
//   double r = 8.0;
  
//   // Get the number of processes in MPI_COMM_WORLD
//   int numProcs;
//   MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

//   // get the rank of this process in MPI_COMM_WORLD
//   int procId;
//   MPI_Comm_rank(MPI_COMM_WORLD, &procId);
  
//   // read in the xyz coordinates
//   int N;
//   double L;
//   std::vector < std::vector < double > > positions = read_xyz_file("data/argon120.xyz", N, L);

//   if (procId == 0){
//     std::cout << "Starting process." << std::endl;
//   }
  
//   double totalCount = 0;
//   std::vector<int> neighbors(N,0);
 
//   bruteForceWithMPI(positions, neighbors, procId, numProcs, N, r);

//   double count = 0;

//   for (int i=0; i < N; i++){
//     count += neighbors[i];
//   }
  
//   MPI_Reduce(&count, &totalCount, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);


//   if (procId==0){
//     std::cout << totalCount/N << std::endl;
//   }

//   // finalise MPI
//   MPI_Finalize();
  
//   return 0;
// }

// void bruteForceWithMPI(std::vector < std::vector < double > > positions, std::vector<int>& neighbors, int procId, int numProcs, int N, double r){
    
//     for (int i = procId - 1; i < N-1; i+= (numProcs-1)){
//       for (int j=i+1; j<N; j++){
//         if (calc_square_distance(positions[i], positions[j])  < (r * r)) {
//           neighbors[i] += 1;
//           neighbors[j] += 1;
//         }
//       }
//     }    
// }



int main(int argc, char **argv)
{
  // initialise MPI
  MPI_Init(&argc, &argv);

  //Choose r
  double r = 8.0;
  
  // Get the number of processes in MPI_COMM_WORLD
  int numProcs;
  MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

  // get the rank of this process in MPI_COMM_WORLD
  int procId;
  MPI_Comm_rank(MPI_COMM_WORLD, &procId);
  
  // read in the xyz coordinates
  int N;
  double L;
  std::string filename = "argon120.xyz";
  std::vector < std::vector < double > > positions = read_xyz_file(filename, N, L);

  //What we add all to
  std::vector <int> totalNeighbors(N,0);
  //Add in a timer 
  double startTime = MPI_Wtime();

  //On master node
  if (procId == 0)
  { 
    std::cout << "Reading in xyz coordinates from " << filename << std::endl;
    std::cout << " Calculating neighbours with " << numProcs << " processes" << std::endl;

    //Initialise the i value and worker (1 because 0 is master)
    int i = 0;
    int worker = 1;

    //While there are spare processes !!!
    while (worker < numProcs)
    {
      //Sending worker an i then update i and worker
      MPI_Send(&i,1,MPI_INT, worker, 0, MPI_COMM_WORLD); 
      worker++;
      i++;
    }
  
    // Run until recived vectors from all i !!!
    while (i < (N))
    {
      MPI_Status stat;
      
      // Initialise Recv from worker meaning done
      int complete;
      MPI_Recv(&complete, 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &stat); //Reciveing an array from procs
      
      //See which worker sent this reply
      worker = stat.MPI_SOURCE; 
      
      //If havent sent out all i, send another task to the worker from which recv
      MPI_Send(&i, 1, MPI_INT, worker, 0 ,MPI_COMM_WORLD);
      i++; 
    }

    //No tasks left so want worker to stop I==N
    int stop = -1;
    //Sends -1 to all workers which breaks while loop in workers
    for (int i = 1; i < numProcs; i++){
      MPI_Send(&stop, 1, MPI_INT, i,0,MPI_COMM_WORLD);
    }
  }

  //Create array for data of worker
  std::vector<int> procNeighbors(N, 0);
  
  //Worker processes
  if (procId != 0)
  {
    // Initialise the i and stat
    int i;
    MPI_Status stat;    

    //Receive the i to work on from master node
    MPI_Recv(&i, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &stat); // Get first task

    //While there is work i!= -1, do the calculation
    while (i != -1)
    {
        //Run bruteForce function, puts data into array
        bruteForce(positions, procNeighbors, N, i, r); 
      
        //Send confirmation and get next i
        int confirm = 1;
        MPI_Send(&confirm, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
      
        MPI_Recv(&i, 1, MPI_INT, 0, 0,MPI_COMM_WORLD, &stat);
    }
  }

  MPI_Reduce(procNeighbors.data(),totalNeighbors.data(), N, MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);

  if (procId == 0){    
    // Print out values
    double average_neighbors;
    int max_neighbors;
    int min_neighbors;
    double timeTaken = MPI_Wtime() - startTime;
    
    getMinMaxAvg(totalNeighbors, N, max_neighbors, min_neighbors, average_neighbors);

    std::cout << "Average neighbors: " << average_neighbors << std::endl;
    std::cout << "Max neighbors: " << max_neighbors << std::endl;
    std::cout << "Min neighbors: " << min_neighbors << std::endl;
    std::cout << "Wall time : " << timeTaken << std::endl;
    std::cout << "Process Over. " << std::endl;
    std::cout << "_____________________ " << std::endl;
  }

  // finalise MPI
  MPI_Finalize();
  
  return 0;
}


void bruteForce(std::vector < std::vector < double > > positions, std::vector<int>& neighbors, int N, int i, double r){
    
    for (int j=i+1; j<N; j++){
      if (calc_square_distance(positions[i], positions[j])  < (r * r)) {
        neighbors[i] += 1;
        neighbors[j] += 1;
      }
    }
}


void bruteForceWithoutMPI(std::vector < std::vector < double > > positions, int N, double L, double r){
  
  std::vector<int> neighbors(N, 0);
  
  for (int i=0; i<N; i++){ //For each atom
    
    for (int j=i+1; j<N; j++){
      if (calc_square_distance(positions[i], positions[j])  < (r * r)) {
        neighbors[i] += 1;
        neighbors[j] += 1;
      }
    }
  }
    
  double average_neighbors;
  int max_neighbors;
  int min_neighbors;
   
  getMinMaxAvg(neighbors, N, max_neighbors, min_neighbors, average_neighbors);

  std::cout << "Average neighbors: " << average_neighbors << std::endl;
  std::cout << "Max neighbors: " << max_neighbors << std::endl;
  std::cout << "Min neighbors: " << min_neighbors << std::endl;

}


void getMinMaxAvg(std::vector<int>& arr , int N, int &max, int &min, double &avg){ //need a const maybe
   max = arr[0], min = arr[0]; int count = 0;
   for(int i = 1; i < N; i++){
      if(max < arr[i])
         max = arr[i];
      if(min > arr[i])
         min = arr[i];
      count += arr[i];
   }
   avg = (double) count / (double) N;
}


double calc_square_distance(std::vector<double> pos_1, std::vector<double> pos_2){
  double x = pos_1[0] - pos_2[0];
  double y = pos_1[1] - pos_2[1];
  double z = pos_1[2] - pos_2[2];

  double d_squared = x*x + y*y + z*z; 
  return d_squared; 
}

// function to read in a list of 3D coordinates from an .xyz file
// input: the name of the file
std::vector < std::vector < double > > read_xyz_file(std::string filename, int& N, double& L){

 

  // open the file
  std::ifstream xyz_file(filename);

  // read in the number of atoms
  xyz_file >> N;
  
  // read in the cell dimension
  xyz_file >> L;
  
  // now read in the positions, ignoring the atomic species
  std::vector < std::vector < double > > positions;
  std::vector < double> pos = {0, 0, 0};
  std::string dummy; 
  for (int i=0;i<N;i++){
    xyz_file >> dummy >> pos[0] >> pos[1] >> pos[2];
    positions.push_back(pos);           
  }
  
  // close the file
  xyz_file.close();
  
  return positions;
}