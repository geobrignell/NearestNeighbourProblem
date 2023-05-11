#include <vector>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <numeric>
#include <algorithm>
#include <stdexcept>
#include <mpi.h>
#include <string>

//Declarations

//Read in
std::vector < std::vector < double > > read_xyz_file(std::string filename, int& N, double& L);

//Cell list Methods
void noMPICellList(int argc, char **argv, std::string filename, double r, double cellWidth);
void MPICellList(int argc, char **argv, std::string filename, double r, double w);

//bruteForce Methods
void noMPIBruteForce(int argc, char **argv, std::string filename, double r);  //BruteForce no MPI
void noMPIBruteForceOptimised(int argc, char **argv, std::string filename, double r);

void MPIBruteForce(int argc, char **argv, std::string filename, double r);
void MPIBruteForceOptimised(int argc, char **argv, std::string filename, double r);
void MPIBruteForceWorkBalanced(int argc, char **argv, std::string filename, double r);
void MPIBruteForceWorkerMaster(int argc, char **argv, std::string filename, double r);

//Calculation Methods
double calc_square_distance(std::vector<double> pos_1, std::vector<double> pos_2);
void getMinMaxAvg(std::vector<int>& arr , int N, int &max, int &min, double &avg);
void bruteForce(std::vector < std::vector < double > > positions, std::vector<int>& neighbors, int N, int i, double r);
void neighbourCheck(std::vector<std::vector<std::vector<std::vector<int>>>> cellList, std::vector<int>& neighbors, std::vector < std::vector < double > > positions, int x1, int y1, int z1, int x2, int y2, int z2, double r);
void selfCheck(std::vector<std::vector<std::vector<std::vector<int>>>> cellList, std::vector<int>& neighbors, std::vector < std::vector < double > > positions, int x1, int y1, int z1, double r);

// =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

//Main 
int main(int argc, char **argv){

  //Change the parameters here:  
  std::string filename = "argon120.xyz"; //File names are argon{num}.xyz, num = 120, 10549, 147023
  double r = 8.0;
  double cellWidth = 9.0;
  
  //Methods: (Uncomment the method you would like to use.)

  noMPICellList(argc, argv, filename, r, cellWidth);
  //MPICellList(argc, argv, filename, r, cellWidth);
  
  //noMPIBruteForce(argc, argv, filename, r);
  //noMPIBruteForceOptimised(argc, argv, filename, r);

  //MPIBruteForce(argc, argv, filename, r);
  //MPIBruteForceOptimised(argc, argv, filename, r);
  //MPIBruteForceWorkBalanced(argc, argv, filename, r);
  //MPIBruteForceWorkerMaster(argc, argv, filename,r );
  
  return 0;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//Cell list Methods

//Cell list using MPI
void MPICellList(int argc, char **argv, std::string filename, double r, double cellWidth){

  // initialise MPI
  MPI_Init(&argc, &argv);

  int procId;
  MPI_Comm_rank(MPI_COMM_WORLD, &procId);

  int numProcs;
  MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

  //Read in File
  int N; double L;
  std::vector < std::vector < double > > positions = read_xyz_file(filename, N, L);

  int numCells = L/cellWidth;

  //Creating the 4D std::vector where first 3 dimensions are of fixed length numCells. Last dimension is of variable size for the atom indexes.
  std::vector<std::vector<std::vector<std::vector<int>>>> cellList(numCells, std::vector<std::vector<std::vector<int>>>(numCells,std::vector<std::vector<int>>(numCells,std::vector<int>() )));

  //Start MPI clock
  double startTime = MPI_Wtime();

  //Putting all atoms into the cell list. Cell list has an array of atom indexs at [ix][iy][iz]
  for(int i=0; i < N; i++){
    int ix = positions[i][0]/cellWidth;
    int iy = positions[i][1]/cellWidth;
    int iz = positions[i][2]/cellWidth;

    //This puts atoms on the far side boundary into the closest cell (e.g. z=18.000 for argon120.xyz)
    if (ix == numCells) ix--;
    if (iy == numCells) iy--;
    if (iz == numCells) iz--;

    cellList[ix][iy][iz].push_back(i);
  }

  //Initialise neighbour count vectors
  std::vector<int> neighbors(N, 0);
  std::vector<int> totalNeighbors(N, 0);

  //For processes not on the main node
  if (procId != 0){
    
    //Looping through all cells
    for (int x= procId-1; x<numCells; x+=(numProcs-1)){
      for (int y= 0; y<numCells; y++){
        for (int z= 0; z<numCells; z++){

          //Loop over all neighbouring cells 
          if (z-1 >= 0){
            if (y-1 >= 0){

              if (x-1 >= 0){
                // then x-1, y-1, z-1 is a neighbour
                //do calculations between cells, add to neighbour list
                neighbourCheck(cellList, neighbors, positions, x, y, z, x-1, y-1, z-1, r);
              }
              if (x+1 < numCells){
                //then x+1, y-1, z-1 is a neighbour
                neighbourCheck(cellList, neighbors, positions, x, y, z, x+1, y-1, z-1, r);
              }
              // x, y-1, z-1 is a neighbour
              neighbourCheck(cellList, neighbors, positions, x, y, z, x, y-1, z-1, r);
            }

            if (y+1 < numCells){

              if (x-1 >= 0){
                // then x-1, y+1, z-1 is a neighbour
                neighbourCheck(cellList, neighbors, positions, x, y, z, x-1, y+1, z-1, r);
              }
              if (x+1 < numCells){
                //then x+1, y+1, z-1 is a neighbour
                neighbourCheck(cellList, neighbors, positions, x, y, z, x+1, y+1, z-1, r);
              }
              // x, y+1, z-1 is a neighbour
              neighbourCheck(cellList, neighbors, positions, x, y, z, x, y+1, z-1, r);
            }

            if (x-1 >= 0){
              // then x-1, y, z-1 is a neighbour
              neighbourCheck(cellList, neighbors, positions, x, y, z, x-1, y, z-1, r);
            }
            if (x+1 < numCells){
              //then x+1, y, z-1 is a neighbour
              neighbourCheck(cellList, neighbors, positions, x, y, z, x+1, y, z-1, r);
            }
            // x, y, z-1 is a neighbour
            neighbourCheck(cellList, neighbors, positions, x, y, z, x, y, z-1, r);
          }
          
          if (y-1 >= 0){

            if (x-1 >= 0){
              // then x-1, y-1, z is a neighbour
              neighbourCheck(cellList, neighbors, positions, x, y, z, x-1, y-1, z, r);
            }
            if (x+1 < numCells){
              //then x+1, y-1, z is a neighbour
              neighbourCheck(cellList, neighbors, positions, x, y, z, x+1, y-1, z, r);
            }
            // x, y-1, z is a neighbour
            neighbourCheck(cellList, neighbors, positions, x, y, z, x, y-1, z, r);
          }

          if (y+1 < numCells){

            if (x-1 >= 0){
              // then x-1, y+1, z is a neighbour
              neighbourCheck(cellList, neighbors, positions, x, y, z, x-1, y+1, z, r);
            }
            if (x+1 < numCells){
              //then x+1, y+1, z is a neighbour
              neighbourCheck(cellList, neighbors, positions, x, y, z, x+1, y+1, z, r);
            }
            // x, y+1, z is a neighbour
            neighbourCheck(cellList, neighbors, positions, x, y, z, x, y+1, z, r);
          }

          if (x-1 >= 0){
            // then x-1, y+1, z is a neighbour
            neighbourCheck(cellList, neighbors, positions, x, y, z, x-1, y, z, r);
          }
          if (x+1 < numCells){
            //then x+1, y+1, z is a neighbour
            neighbourCheck(cellList, neighbors, positions, x, y, z, x+1, y, z, r);
          }
            
          if ( z+1 < numCells){
            if (y-1 >= 0){

              if (x-1 >= 0){
                // then x-1, y-1, z+1 is a neighbour
                neighbourCheck(cellList, neighbors, positions, x, y, z, x-1, y-1, z+1, r);
              }
              if (x+1 < numCells){
                //then x+1, y-1, z+1 is a neighbour
                neighbourCheck(cellList, neighbors, positions, x, y, z, x+1, y-1, z+1, r);
              }
              // x, y-1, z+1 is a neighbour
              neighbourCheck(cellList, neighbors, positions, x, y, z, x, y-1, z+1, r);
            }

            if (y+1 < numCells){

              if (x-1 >= 0){
                // then x-1, y+1, z+1 is a neighbour
                neighbourCheck(cellList, neighbors, positions, x, y, z, x-1, y+1, z+1, r);
              }
              if (x+1 < numCells){
                //then x+1, y+1, z+1 is a neighbour
                neighbourCheck(cellList, neighbors, positions, x, y, z, x+1, y+1, z+1, r);
              }
              // x, y+1, z+1 is a neighbour
              neighbourCheck(cellList, neighbors, positions, x, y, z, x, y+1, z+1, r);
            }

            if (x-1 >= 0){
              // then x-1, y, z+1 is a neighbour
              neighbourCheck(cellList, neighbors, positions, x, y, z, x-1, y, z+1, r);
            }
            if (x+1 < numCells){
              //then x+1, y, z+1 is a neighbour
              neighbourCheck(cellList, neighbors, positions, x, y, z, x+1, y, z+1, r);
            }
            // x, y, z+1 is a neighbour
            neighbourCheck(cellList, neighbors, positions, x, y, z, x, y, z+1, r);
          }

          // Own cell
          selfCheck(cellList, neighbors, positions, x, y, z, r);
        }
      }
    }
  }

  //Reduce to main node
  MPI_Reduce(neighbors.data(),totalNeighbors.data(), N, MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);

  //Print Outputs
  if (procId == 0){
    std::cout << "Using the MPICellList method. " << std::endl;
    std::cout << "Reading in xyz coordinates from " << filename << std::endl;
    std::cout << "Calculating neighbours with " << numProcs << " processes." << std::endl;
    std::cout << "Number of cells : " << numCells << ", Length L : " << L << ", Width : " << cellWidth << "." << std::endl;

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
  
  MPI_Finalize();
}


//Cell list without MPI (1 process)
void noMPICellList(int argc, char **argv, std::string filename, double r, double cellWidth){

  // initialise MPI
  MPI_Init(&argc, &argv);

  int procId;
  MPI_Comm_rank(MPI_COMM_WORLD, &procId);

  int numProcs;
  MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

  //As we only want one process to run. Does everything on a single node.
  if (procId == 0){ 

    //Reads in file
    int N; double L;
    std::vector < std::vector < double > > positions = read_xyz_file(filename, N, L);
    double startTime = MPI_Wtime();

    int numCells = L/cellWidth;
    std::cout << "Number of cells : " << numCells << std::endl;

    //Creating the 4D std::vector where the 3 dimensions are fixed length of numCells and the last one is variable for the atom index to be put into
    std::vector<std::vector<std::vector<std::vector<int>>>> cellList(numCells, std::vector<std::vector<std::vector<int>>>(numCells,std::vector<std::vector<int>>(numCells,std::vector<int>() )));

    //Putting all atoms into the cell list. Cell list has an array of atom indexs
    for(int i=0; i < N; i++){
      int ix = positions[i][0]/cellWidth;
      int iy = positions[i][1]/cellWidth;
      int iz = positions[i][2]/cellWidth;

      if (ix == numCells) ix--;
      if (iy == numCells) iy--;
      if (iz == numCells) iz--;

      cellList[ix][iy][iz].push_back(i);
    }

    //std::vector to store the number of neighbours for i
    std::vector<int> neighbors(N, 0);

    //std::cout << "Looping through cells" << std::endl;
    //Looping through all cells
    for (int x=0; x<numCells; x++){
      for (int y=0; y<numCells; y++){
        for (int z=0; z<numCells; z++){

          
          //std::cout << "Cell " << x << ", " << y << ", " << z << std::endl;
          //Loop over all neighbouring cells 
          if (z-1 >= 0){
            if (y-1 >= 0){

              if (x-1 >= 0){
                // then x-1, y-1, z-1 is a neighbour
                //do calculations between cells, add to neighbour list
                neighbourCheck(cellList, neighbors, positions, x, y, z, x-1, y-1, z-1, r);
              }

              if (x+1 < numCells){
                //then x+1, y-1, z-1 is a neighbour
                neighbourCheck(cellList, neighbors, positions, x, y, z, x+1, y-1, z-1, r);
              }
              // x, y-1, z-1 is a neighbour
              neighbourCheck(cellList, neighbors, positions, x, y, z, x, y-1, z-1, r);
            }

            if (y+1 < numCells){

              if (x-1 >= 0){
                // then x-1, y+1, z-1 is a neighbour
                neighbourCheck(cellList, neighbors, positions, x, y, z, x-1, y+1, z-1, r);
              }
              if (x+1 < numCells){
                //then x+1, y+1, z-1 is a neighbour
                neighbourCheck(cellList, neighbors, positions, x, y, z, x+1, y+1, z-1, r);
              }
              // x, y+1, z-1 is a neighbour
              neighbourCheck(cellList, neighbors, positions, x, y, z, x, y+1, z-1, r);
            }

            if (x-1 >= 0){
              // then x-1, y, z-1 is a neighbour
              neighbourCheck(cellList, neighbors, positions, x, y, z, x-1, y, z-1, r);
            }
            if (x+1 < numCells){
              //then x+1, y, z-1 is a neighbour
              neighbourCheck(cellList, neighbors, positions, x, y, z, x+1, y, z-1, r);
            }
            // x, y, z-1 is a neighbour
            neighbourCheck(cellList, neighbors, positions, x, y, z, x, y, z-1, r);
            
          }
          
          if (y-1 >= 0){

            if (x-1 >= 0){
              // then x-1, y-1, z is a neighbour
              neighbourCheck(cellList, neighbors, positions, x, y, z, x-1, y-1, z, r);
            }
            if (x+1 < numCells){
              //then x+1, y-1, z is a neighbour
              neighbourCheck(cellList, neighbors, positions, x, y, z, x+1, y-1, z, r);
            }
            // x, y-1, z is a neighbour
            neighbourCheck(cellList, neighbors, positions, x, y, z, x, y-1, z, r);
          }

          if (y+1 < numCells){

            if (x-1 >= 0){
              // then x-1, y+1, z is a neighbour
              neighbourCheck(cellList, neighbors, positions, x, y, z, x-1, y+1, z, r);
            }
            if (x+1 < numCells){
              //then x+1, y+1, z is a neighbour
              neighbourCheck(cellList, neighbors, positions, x, y, z, x+1, y+1, z, r);
            }
            // x, y+1, z is a neighbour
            neighbourCheck(cellList, neighbors, positions, x, y, z, x, y+1, z, r);
          }

          if (x-1 >= 0){
            // then x-1, y+1, z is a neighbour
            neighbourCheck(cellList, neighbors, positions, x, y, z, x-1, y, z, r);
          }
          if (x+1 < numCells){
            //then x+1, y+1, z is a neighbour
            neighbourCheck(cellList, neighbors, positions, x, y, z, x+1, y, z, r);
          }
            

          if ( z+1 < numCells){
            if (y-1 >= 0){

              if (x-1 >= 0){
                // then x-1, y-1, z+1 is a neighbour
                neighbourCheck(cellList, neighbors, positions, x, y, z, x-1, y-1, z+1, r);
              }
              if (x+1 < numCells){
                //then x+1, y-1, z+1 is a neighbour
                neighbourCheck(cellList, neighbors, positions, x, y, z, x+1, y-1, z+1, r);
              }
              // x, y-1, z+1 is a neighbour
              neighbourCheck(cellList, neighbors, positions, x, y, z, x, y-1, z+1, r);
            }

            if (y+1 < numCells){

              if (x-1 >= 0){
                // then x-1, y+1, z+1 is a neighbour
                neighbourCheck(cellList, neighbors, positions, x, y, z, x-1, y+1, z+1, r);
              }
              if (x+1 < numCells){
                //then x+1, y+1, z+1 is a neighbour
                neighbourCheck(cellList, neighbors, positions, x, y, z, x+1, y+1, z+1, r);
              }
              // x, y+1, z+1 is a neighbour
              neighbourCheck(cellList, neighbors, positions, x, y, z, x, y+1, z+1, r);
            }

          
            if (x-1 >= 0){
              // then x-1, y, z+1 is a neighbour
              neighbourCheck(cellList, neighbors, positions, x, y, z, x-1, y, z+1, r);
            }
            if (x+1 < numCells){
              //then x+1, y, z+1 is a neighbour
              neighbourCheck(cellList, neighbors, positions, x, y, z, x+1, y, z+1, r);
            }
            // x, y, z+1 is a neighbour
            neighbourCheck(cellList, neighbors, positions, x, y, z, x, y, z+1, r);
        
          }

          // Own cell
          selfCheck(cellList, neighbors, positions, x, y, z, r);

        }
      }
    }

    //Printing outputs
    std::cout << "Reading in xyz coordinates from " << filename << std::endl;
    std::cout << " Calculating neighbours with " << numProcs << " processes" << std::endl;

    double average_neighbors;
    int max_neighbors;
    int min_neighbors;
    double timeTaken = MPI_Wtime() - startTime;

      
    getMinMaxAvg(neighbors, N, max_neighbors, min_neighbors, average_neighbors);

    std::cout << "Using the noMPICellList method. " << std::endl;
    std::cout << "Average neighbors: " << average_neighbors << std::endl;
    std::cout << "Max neighbors: " << max_neighbors << std::endl;
    std::cout << "Min neighbors: " << min_neighbors << std::endl;
    std::cout << "Wall time : " << timeTaken << std::endl;
    std::cout << "Process Over. " << std::endl;
    std::cout << "_____________________ " << std::endl;
  }
  
  MPI_Finalize();
}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


//Brute force methods
// Basic brute force split by processes (j=0)
void MPIBruteForce(int argc, char **argv, std::string filename, double r){
  
  // initialise MPI
  MPI_Init(&argc, &argv);
  
  // Get the number of processes in MPI_COMM_WORLD
  int numProcs;
  MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

  // get the rank of this process in MPI_COMM_WORLD
  int procId;
  MPI_Comm_rank(MPI_COMM_WORLD, &procId);
  
  // Read in the xyz coordinates
  int N;
  double L;
  std::vector < std::vector < double > > positions = read_xyz_file(filename, N, L);
  double startTime = MPI_Wtime();

  //Initialise neighbour count vectors
  std::vector <int> totalNeighbors(N,0);
  std::vector <int> localNeighbors(N,0);

  //Loop through all atoms in file with every other atom in file and add to neighbour if within r.
  for (int i= procId; i<N; i+= numProcs){
      for (int j=0; j<N; j++){
        if (i!=j){
          if (calc_square_distance(positions[i], positions[j])  < (r * r)) {
              localNeighbors[i] += 1;
          }
        }
      }
  }

  //Reduce all to root node.
  MPI_Reduce(localNeighbors.data(),totalNeighbors.data(), N, MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
  
  // Print out values
  if (procId == 0)
  {
  std::cout << "Reading in xyz coordinates from " << filename << std::endl;
  std::cout << " Calculating neighbours with " << numProcs << " processes" << std::endl;
  
  double average_neighbors;
  int max_neighbors;
  int min_neighbors;
  double timeTaken = MPI_Wtime() - startTime;
  
  getMinMaxAvg(totalNeighbors, N, max_neighbors, min_neighbors, average_neighbors);
  
  std::cout << "Using the MPIBruteForce method. " << std::endl;
  std::cout << "Average neighbors: " << average_neighbors << std::endl;
  std::cout << "Max neighbors: " << max_neighbors << std::endl;
  std::cout << "Min neighbors: " << min_neighbors << std::endl;
  std::cout << "Wall time : " << timeTaken << std::endl;
  std::cout << "Process Over. " << std::endl;
  std::cout << "_____________________ " << std::endl;
  }

  MPI_Finalize();
}


// Basic brute force split by processes but optimised using j = i+1
void MPIBruteForceOptimised(int argc, char **argv, std::string filename, double r){
  
  // initialise MPI
  MPI_Init(&argc, &argv);
  
  // Get the number of processes in MPI_COMM_WORLD
  int numProcs;
  MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

  // get the rank of this process in MPI_COMM_WORLD
  int procId;
  MPI_Comm_rank(MPI_COMM_WORLD, &procId);
  
  // Read in the xyz coordinates
  int N;
  double L;
  std::vector < std::vector < double > > positions = read_xyz_file(filename, N, L);
  double startTime = MPI_Wtime();

  //Initialise neighbour count vectors
  std::vector <int> totalNeighbors(N,0);
  std::vector <int> localNeighbors(N,0);
  
  //Loop through all atoms, with i+1 atoms to reduce load.
  for (int i= procId; i<N; i+= numProcs){
      for (int j=i+1; j<N; j++){
          if (calc_square_distance(positions[i], positions[j])  < (r * r)) {
              localNeighbors[i] += 1;
              localNeighbors[j] += 1;
          }
      }
  }
  
  //Reduce to root node
  MPI_Reduce(localNeighbors.data(),totalNeighbors.data(), N, MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);

  // Printing out values
  if (procId == 0)
  {
  std::cout << "Reading in xyz coordinates from " << filename << std::endl;
  std::cout << " Calculating neighbours with " << numProcs << " processes" << std::endl;

  double average_neighbors;
  int max_neighbors;
  int min_neighbors;
  double timeTaken = MPI_Wtime() - startTime;
  
  getMinMaxAvg(totalNeighbors, N, max_neighbors, min_neighbors, average_neighbors);
  
  std::cout << "Using the BruteForceOptimised method. " << std::endl;
  std::cout << "Average neighbors: " << average_neighbors << std::endl;
  std::cout << "Max neighbors: " << max_neighbors << std::endl;
  std::cout << "Min neighbors: " << min_neighbors << std::endl;
  std::cout << "Wall time : " << timeTaken << std::endl;
  std::cout << "Process Over. " << std::endl;
  std::cout << "_____________________ " << std::endl;
  }

  // finalise MPI
  MPI_Finalize();
}

void MPIBruteForceWorkBalanced(int argc, char **argv, std::string filename, double r){
  
  // initialise MPI
  MPI_Init(&argc, &argv);
  
  // Get the number of processes in MPI_COMM_WORLD
  int numProcs;
  MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

  // get the rank of this process in MPI_COMM_WORLD
  int procId;
  MPI_Comm_rank(MPI_COMM_WORLD, &procId);
  
  // read in the xyz coordinates
  int N;
  double L;
  std::vector < std::vector < double > > positions = read_xyz_file(filename, N, L);

  double startTime = MPI_Wtime();
  
  //Initalise neighbour count vectors
  std::vector <int> totalNeighbors(N,0);
  std::vector <int> localNeighbors(N,0);

  //Calculates the amount of work and work per process
  long long int totJobs = 0.5*N * (N-1);
  long long int jobsPerProc = floor(totJobs/numProcs);

  //Initialises iLimits which will store the starting i values for all the processes
  std::vector <int> iLimits;

  //Function that appends to iLimits when the amount of work is above the work per process value
  long long int count = 0;
  for (int i = N-1; i>0; --i){
    count += i;
    if(count >= jobsPerProc){
      iLimits.push_back(i);
      count = 0;
    }
  }

  // Boundary variables for each process
  int i0; int i1;

  //Root node does from i=0
  if (procId ==0){
    i0 = 0;
    i1 = N - iLimits[procId];
  }

  //If not the root node, gets its boundary values from iLimits vector
  else {
    i0 = N - iLimits[procId -1];
    i1 = N - iLimits[procId];
  }

  //Does the brute force between its unique limits
  for (int i = i0; i < i1; i++){
    for (int j=i+1; j<N; j++){ 
      if (calc_square_distance(positions[i], positions[j])  < (r * r)) {
        localNeighbors[i] += 1;
        localNeighbors[j] += 1;
      }
    }
  }

  //Reducing to root node.
  MPI_Reduce(localNeighbors.data(),totalNeighbors.data(), N, MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
  
  //Printing outputs
  if (procId == 0)
  {
  std::cout << "Reading in xyz coordinates from " << filename << std::endl;
  std::cout << " Calculating neighbours with " << numProcs << " processes" << std::endl;
  double average_neighbors;
  int max_neighbors;
  int min_neighbors;
  double timeTaken = MPI_Wtime() - startTime;
  
  getMinMaxAvg(totalNeighbors, N, max_neighbors, min_neighbors, average_neighbors);
  
  std::cout << "Using the MPIBruteForceWorkBalanced method. " << std::endl;
  std::cout << "Average neighbors: " << average_neighbors << std::endl;
  std::cout << "Max neighbors: " << max_neighbors << std::endl;
  std::cout << "Min neighbors: " << min_neighbors << std::endl;
  std::cout << "Wall time : " << timeTaken << std::endl;
  std::cout << "Process Over. " << std::endl;
  std::cout << "_____________________ " << std::endl;
  }

  MPI_Finalize();
}

// Brute force with workload managed using worker-manager system
void MPIBruteForceWorkerMaster(int argc, char **argv, std::string filename, double r)
{
  // initialise MPI
  MPI_Init(&argc, &argv);
  
  // Get the number of processes in MPI_COMM_WORLD
  int numProcs;
  MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

  // get the rank of this process in MPI_COMM_WORLD
  int procId;
  MPI_Comm_rank(MPI_COMM_WORLD, &procId);
  
  // read in the xyz coordinates
  int N;
  double L;
  std::vector < std::vector < double > > positions = read_xyz_file(filename, N, L);
  double startTime = MPI_Wtime();
  
  //Vector in master node to be reduced to
  std::vector <int> totalNeighbors(N,0);
  
  //On master node
  if (procId == 0)
  { 
    //Initialise the i value and worker (1 because 0 is master)
    int i = 0;
    int worker = 1;

    //While there are spare processes
    while (worker < numProcs)
    {
      //Sending worker an i then update i and worker
      MPI_Send(&i,1,MPI_INT, worker, 0, MPI_COMM_WORLD); 
      worker++;
      i++;
    }
  
    // Run until all i values have been distributed
    while (i < (N))
    {
      MPI_Status stat;
      
      // Recv from worker meaning that they're done
      int complete;
      MPI_Recv(&complete, 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &stat); //Reciveing an array from procs
      
      //See which worker sent this reply
      worker = stat.MPI_SOURCE; 
      
      //If havent sent out all i, send another task to the worker from which recveived
      MPI_Send(&i, 1, MPI_INT, worker, 0 ,MPI_COMM_WORLD);
      i++; 
    }

    //No tasks left so want worker to stop i==N
    int stop = -1;
    //Sends -1 to all workers which breaks while loop in workers
    for (int i = 1; i < numProcs; i++){
      MPI_Send(&stop, 1, MPI_INT, i,0,MPI_COMM_WORLD);
    }
  }

  //Create local vector for worker
  std::vector<int> procNeighbors(N, 0);
  
  //Worker processes
  if (procId != 0)
  {
    // Initialise the i and stat
    int i;
    MPI_Status stat;    

    //Receive the i to work on from master node
    MPI_Recv(&i, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &stat); // Get first task

    //While there is work, i.e. i!= -1, do the calculation
    while (i != -1)
    {
        //Run bruteForce function, puts data into local vector
        bruteForce(positions, procNeighbors, N, i, r); 
      
        //Send confirmation and get next i
        int confirm = 1;
        MPI_Send(&confirm, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);

        MPI_Recv(&i, 1, MPI_INT, 0, 0,MPI_COMM_WORLD, &stat);
    }
  }
  //Reduce to master node
  MPI_Reduce(procNeighbors.data(),totalNeighbors.data(), N, MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);

  //Print outputs
  if (procId == 0){    
    std::cout << "Reading in xyz coordinates from " << filename << std::endl;
    std::cout << " Calculating neighbours with " << numProcs << " processes" << std::endl;
    
    double average_neighbors;
    int max_neighbors;
    int min_neighbors;
    double timeTaken = MPI_Wtime() - startTime;
    
    getMinMaxAvg(totalNeighbors, N, max_neighbors, min_neighbors, average_neighbors);

    std::cout << "Using the MPIBruteForceWorkloadManaged method. " << std::endl;
    std::cout << "Average neighbors: " << average_neighbors << std::endl;
    std::cout << "Max neighbors: " << max_neighbors << std::endl;
    std::cout << "Min neighbors: " << min_neighbors << std::endl;
    std::cout << "Wall time : " << timeTaken << std::endl;
    std::cout << "Process Over. " << std::endl;
    std::cout << "_____________________ " << std::endl;
  }

  MPI_Finalize();  
}

//BruteForce without MPI (j=0)
void noMPIBruteForce(int argc, char **argv, std::string filename, double r){
  
  // initialise MPI
  MPI_Init(&argc, &argv);

  int procId;
  MPI_Comm_rank(MPI_COMM_WORLD, &procId);

  int numProcs;
  MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

  //Only want one process
  if (procId == 0){

    //Read in the xyz coordinates
    int N; double L;
    std::vector < std::vector < double > > positions = read_xyz_file(filename, N, L);
    double startTime = MPI_Wtime();

    //Create a vector for neighbour data
    std::vector<int> neighbors(N, 0);
    
    //Loop through all atoms with every other atom and add values to vector
    for (int i=0; i<N; i++){
      
      for (int j=0; j<N; j++){
        if (i!=j){
          if (calc_square_distance(positions[i], positions[j])  < (r * r)) {
            neighbors[i] += 1;
          }
        }
      }
    }
    
    //Print outputs
    std::cout << "Reading in xyz coordinates from " << filename << std::endl;
    std::cout << " Calculating neighbours with " << numProcs << " processes" << std::endl;
    double average_neighbors;
    int max_neighbors;
    int min_neighbors;
    double timeTaken = MPI_Wtime() - startTime;
    
    getMinMaxAvg(neighbors, N, max_neighbors, min_neighbors, average_neighbors);

    std::cout << "Using the noMPIBruteForce method. " << std::endl;
    std::cout << "Average neighbors: " << average_neighbors << std::endl;
    std::cout << "Max neighbors: " << max_neighbors << std::endl;
    std::cout << "Min neighbors: " << min_neighbors << std::endl;
    std::cout << "Wall time : " << timeTaken << std::endl;
    std::cout << "Process Over. " << std::endl;
    std::cout << "_____________________ " << std::endl;
  }

  MPI_Finalize();
}

//Function : BruteForce without MPI j=i+1
void noMPIBruteForceOptimised(int argc, char **argv, std::string filename, double r){
  
  // initialise MPI
  MPI_Init(&argc, &argv);

  int procId;
  MPI_Comm_rank(MPI_COMM_WORLD, &procId);

  int numProcs;
  MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

  //Only want one process
  if (procId == 0){ 

    // Read in the xyz coordinates
    int N; double L;
    std::vector < std::vector < double > > positions = read_xyz_file(filename, N, L);
    double startTime = MPI_Wtime();

    // Vector to store neighbour data
    std::vector<int> neighbors(N, 0);
    
    //Loop over every atom with all the other atoms that havent been calculated yet.
    for (int i=0; i<N; i++){ 
      
      for (int j=i+1; j<N; j++){
        if (calc_square_distance(positions[i], positions[j])  < (r * r)) {
          neighbors[i] += 1;
          neighbors[j] += 1;
        }
      }
    }
    
    //Print outputs
    double average_neighbors;
    int max_neighbors;
    int min_neighbors;
    double timeTaken = MPI_Wtime() - startTime;

    std::cout << "Reading in xyz coordinates from " << filename << std::endl;
    std::cout << " Calculating neighbours with " << numProcs << " processes" << std::endl;
    
    getMinMaxAvg(neighbors, N, max_neighbors, min_neighbors, average_neighbors);

    std::cout << "Using the noMPIBruteForceOptimised method. " << std::endl;
    std::cout << "Average neighbors: " << average_neighbors << std::endl;
    std::cout << "Max neighbors: " << max_neighbors << std::endl;
    std::cout << "Min neighbors: " << min_neighbors << std::endl;
    std::cout << "Wall time : " << timeTaken << std::endl;
    std::cout << "Process Over. " << std::endl;
    std::cout << "_____________________ " << std::endl;
  }
  MPI_Finalize();

}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//Calculation Methods

//Method which loops over a cells indexes with a neighbouring cells indexes and performs the brute force method to add to a neighbour count vector
void neighbourCheck(std::vector<std::vector<std::vector<std::vector<int>>>> cellList, std::vector<int>& neighbors, std::vector < std::vector < double > > positions, int x1, int y1, int z1, int x2, int y2, int z2, double r){
  //For all the indexs stored in cell 1, see if in range with all indexs in cell 2
  // if in range, add neighbor count to cell 1 index
  for (int i=0; i<(cellList[x1][y1][z1].size()); i++){
        for (int j=0; j<(cellList[x2][y2][z2].size()); j++)
        {
          if (calc_square_distance(positions[cellList[x1][y1][z1][i]], positions[cellList[x2][y2][z2][j]])  < (r * r))
          {
            neighbors[cellList[x1][y1][z1][i]] += 1;
          }
        }
  }
}

// Brute forces all the atoms within a cell with all the others.
void selfCheck(std::vector<std::vector<std::vector<std::vector<int>>>> cellList, std::vector<int>& neighbors, std::vector < std::vector < double > > positions, int x1, int y1, int z1, double r){
  // if in range, add a neighbor count to cell 1 index
  int size = cellList[x1][y1][z1].size();
  for (int i=0; i<size; i++){
    for (int j=0; j<size; j++)
    {
      // Do not want to compare between the same index
      if (i != j){
        if (calc_square_distance(positions[cellList[x1][y1][z1][i]], positions[cellList[x1][y1][z1][j]])  < (r * r))
        {
          neighbors[cellList[x1][y1][z1][i]] += 1;
        }
      }
    }
  }
}


//Optimised brute force method called in the complicated methods
void bruteForce(std::vector < std::vector < double > > positions, std::vector<int>& neighbors, int N, int i, double r){
    
    for (int j=i+1; j<N; j++){
      if (calc_square_distance(positions[i], positions[j])  < (r * r)) {
        neighbors[i] += 1;
        neighbors[j] += 1;
      }
    }
}


// Changes the referenced values to the max min and avg values of a vector
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


// Returns the d**2 value between two positions
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