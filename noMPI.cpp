#include <vector>
#include <stdlib.h>
#include <fstream>

#include <iostream>
#include<numeric>
#include <map>
#include <algorithm>
#include <set>
#include<stdexcept>
//#include <mpi.h>


double calc_square_distance(std::vector<double> pos_1, std::vector<double> pos_2);
std::vector < std::vector < double > > read_xyz_file(std::string filename, int& N, double& L);
void getMinMaxAvg(std::vector<int>& arr , int N, int &max, int &min, double &avg);
void niceFunction(std::vector<std::vector<std::vector<std::vector<int>>>> cellList, std::vector<int>& neighbors, std::vector < std::vector < double > > positions, int x1, int y1, int z1, int x2, int y2, int z2, double r);
void ownNiceFunction(std::vector<std::vector<std::vector<std::vector<int>>>> cellList, std::vector<int>& neighbors, std::vector < std::vector < double > > positions, int x1, int y1, int z1, double r);


int main(int argc, char **argv)
{
  double r = 8.0;
  int N; double L;
  std::vector < std::vector < double > > positions = read_xyz_file("data/argon120.xyz", N, L);

  double cellWidth = 18.0; //Will want to vary this
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

    //std::cout << i << " " << ix << " " << iy << " " << iz << std::endl;
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
              niceFunction(cellList, neighbors, positions, x, y, z, x-1, y-1, z-1, r);
            }

            if (x+1 < numCells){
              //then x+1, y-1, z-1 is a neighbour
              niceFunction(cellList, neighbors, positions, x, y, z, x+1, y-1, z-1, r);
            }
            // x, y-1, z-1 is a neighbour
            niceFunction(cellList, neighbors, positions, x, y, z, x, y-1, z-1, r);
          }

          if (y+1 < numCells){

            if (x-1 >= 0){
              // then x-1, y+1, z-1 is a neighbour
              niceFunction(cellList, neighbors, positions, x, y, z, x-1, y+1, z-1, r);
            }
            if (x+1 < numCells){
              //then x+1, y+1, z-1 is a neighbour
              niceFunction(cellList, neighbors, positions, x, y, z, x+1, y+1, z-1, r);
            }
            // x, y+1, z-1 is a neighbour
            niceFunction(cellList, neighbors, positions, x, y, z, x, y+1, z-1, r);
          }

          if (x-1 >= 0){
            // then x-1, y, z-1 is a neighbour
            niceFunction(cellList, neighbors, positions, x, y, z, x-1, y, z-1, r);
          }
          if (x+1 < numCells){
            //then x+1, y, z-1 is a neighbour
            niceFunction(cellList, neighbors, positions, x, y, z, x+1, y, z-1, r);
          }
          // x, y, z-1 is a neighbour
          niceFunction(cellList, neighbors, positions, x, y, z, x, y, z-1, r);
          
        }
        
        if (y-1 >= 0){

          if (x-1 >= 0){
            // then x-1, y-1, z is a neighbour
            niceFunction(cellList, neighbors, positions, x, y, z, x-1, y-1, z, r);
          }
          if (x+1 < numCells){
            //then x+1, y-1, z is a neighbour
            niceFunction(cellList, neighbors, positions, x, y, z, x+1, y-1, z, r);
          }
          // x, y-1, z is a neighbour
          niceFunction(cellList, neighbors, positions, x, y, z, x, y-1, z, r);
        }

        if (y+1 < numCells){

          if (x-1 >= 0){
            // then x-1, y+1, z is a neighbour
            niceFunction(cellList, neighbors, positions, x, y, z, x-1, y+1, z, r);
          }
          if (x+1 < numCells){
            //then x+1, y+1, z is a neighbour
            niceFunction(cellList, neighbors, positions, x, y, z, x+1, y+1, z, r);
          }
          // x, y+1, z is a neighbour
          niceFunction(cellList, neighbors, positions, x, y, z, x, y+1, z, r);
        }

        if (x-1 >= 0){
          // then x-1, y+1, z is a neighbour
          niceFunction(cellList, neighbors, positions, x, y, z, x-1, y, z, r);
        }
        if (x+1 < numCells){
          //then x+1, y+1, z is a neighbour
          niceFunction(cellList, neighbors, positions, x, y, z, x+1, y, z, r);
        }
          

        if ( z+1 < numCells){
          if (y-1 >= 0){

            if (x-1 >= 0){
              // then x-1, y-1, z+1 is a neighbour
              niceFunction(cellList, neighbors, positions, x, y, z, x-1, y-1, z+1, r);
            }
            if (x+1 < numCells){
              //then x+1, y-1, z+1 is a neighbour
              niceFunction(cellList, neighbors, positions, x, y, z, x+1, y-1, z+1, r);
            }
            // x, y-1, z+1 is a neighbour
            niceFunction(cellList, neighbors, positions, x, y, z, x, y-1, z+1, r);
          }

          if (y+1 < numCells){

            if (x-1 >= 0){
              // then x-1, y+1, z+1 is a neighbour
              niceFunction(cellList, neighbors, positions, x, y, z, x-1, y+1, z+1, r);
            }
            if (x+1 < numCells){
              //then x+1, y+1, z+1 is a neighbour
              niceFunction(cellList, neighbors, positions, x, y, z, x+1, y+1, z+1, r);
            }
            // x, y+1, z+1 is a neighbour
            niceFunction(cellList, neighbors, positions, x, y, z, x, y+1, z+1, r);
          }

        
          if (x-1 >= 0){
            // then x-1, y, z+1 is a neighbour
            niceFunction(cellList, neighbors, positions, x, y, z, x-1, y, z+1, r);
          }
          if (x+1 < numCells){
            //then x+1, y, z+1 is a neighbour
            niceFunction(cellList, neighbors, positions, x, y, z, x+1, y, z+1, r);
          }
          // x, y, z+1 is a neighbour
          niceFunction(cellList, neighbors, positions, x, y, z, x, y, z+1, r);
      
        }

        // Own cell
        ownNiceFunction(cellList, neighbors, positions, x, y, z, r);

      }
    }
  }

  for (int i = 0; i < neighbors.size(); i++){
    std::cout << neighbors[i] << std::endl;
  }

  double average_neighbors;
  int max_neighbors;
  int min_neighbors;
   
  getMinMaxAvg(neighbors, N, max_neighbors, min_neighbors, average_neighbors);

  std::cout << "Average neighbors: " << average_neighbors << std::endl;
  std::cout << "Max neighbors: " << max_neighbors << std::endl;
  std::cout << "Min neighbors: " << min_neighbors << std::endl;

  return 0;
}


void ownNiceFunction(std::vector<std::vector<std::vector<std::vector<int>>>> cellList, std::vector<int>& neighbors, std::vector < std::vector < double > > positions, int x1, int y1, int z1, double r){
  //For all the indexs stored in cell 1, see if in range with all indexs in cell 2
  // if in range, add neighbor count to cell 1 index
  int size = cellList[x1][y1][z1].size();
  for (int i=0; i<size; i++){
      for (int j=0; j<size; j++)
      {
        
        if (i != j){
          if (calc_square_distance(positions[cellList[x1][y1][z1][i]], positions[cellList[x1][y1][z1][j]])  < (r * r))
          {
            neighbors[cellList[x1][y1][z1][i]] += 1;
            //std::cout << "Added to " << cellList[x1][y1][z1][i] << std::endl;
          }
        }
      }
  }
}


void niceFunction(std::vector<std::vector<std::vector<std::vector<int>>>> cellList, std::vector<int>& neighbors, std::vector < std::vector < double > > positions, int x1, int y1, int z1, int x2, int y2, int z2, double r){
  //For all the indexs stored in cell 1, see if in range with all indexs in cell 2
  // if in range, add neighbor count to cell 1 index
  for (int i=0; i<(cellList[x1][y1][z1].size()); i++){
        for (int j=0; j<(cellList[x2][y2][z2].size()); j++)
        {
          if (calc_square_distance(positions[cellList[x1][y1][z1][i]], positions[cellList[x2][y2][z2][j]])  < (r * r))
          {
            neighbors[cellList[x1][y1][z1][i]] += 1;
            //std::cout << "Added to " << cellList[x1][y1][z1][i] << std::endl;
          }

        }
  }
}


void bruteForceWithoutMPI(std::vector < std::vector < double > > positions, int N, double r, std::vector<std::vector<std::vector<std::vector<int>>>> cellList){
  
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


// =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
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

  std::cout << "Reading in xyz coordinates from " << filename << std::endl;

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

