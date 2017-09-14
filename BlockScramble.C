#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/time.h>
#include <time.h>
#include <assert.h>
#include <omp.h>



#include <qmp.h>
#include <mpi.h>
#include <vector>


//#include <Grid/Grid.h>
#include "BlockGeom.h"
#include "Scramble.h"

/**
 * Get current time in milli seconds.
 */
static double
dclock (void)
{
  struct timeval tv;

  gettimeofday (&tv, 0);

  return tv.tv_sec*1000.0 + tv.tv_usec/1000.0;
}

#define MEM_SIZE 81920 
typedef int64_t DATA;

//#define PRINT QMP_printf
#define PRINT printf

int init_QMP(int *argc, char*** argv,
    std::vector<int> &GlobalDim,
    std::vector<int> &GlobalPos
	){

  QMP_status_t status, err;
  QMP_thread_level_t req, prv;

  req = QMP_THREAD_MULTIPLE;
  status = QMP_init_msg_passing (argc, argv, req, &prv);
  if (*argc < 3) {
    if (QMP_is_primary_node())
      fprintf (stderr, "Usage: %s numloops msgsize [-v]\n", argv[0]);
    exit (1);
  }

  if (status != QMP_SUCCESS) {
    QMP_error ("QMP_init failed: %s\n", QMP_error_string(status));
    QMP_abort(1);
  }

    int peNum = QMP_get_number_of_nodes();
    int NDIM = QMP_get_allocated_number_of_dimensions();
    const int *peGrid_t = QMP_get_allocated_dimensions();
    const int *pePos_t = QMP_get_allocated_coordinates();
  if (!QMP_get_node_number())
   QMP_printf("QMP threading %d requested, %d provided\n",req,prv);
//    assert(NDIM==4);

//  status = QMP_declare_logical_topology ((int *)GlobalDim.data(), NDIM);
  status = QMP_declare_logical_topology (peGrid_t, NDIM);
  if (status != QMP_SUCCESS)
    PRINT ("Cannot declare logical grid\n");

GlobalDim.resize(NDIM);
GlobalPos.resize(NDIM);

    int *SplitGrid  = new int[NDIM];
    int *SplitPos  = new int[NDIM];


    for(int i =0;i<NDIM;i++){
        GlobalDim[i] = peGrid_t[i];
        GlobalPos[i] = pePos_t[i];
    }

	return NDIM;
}

#if 0
int init_Grid(int *argc, char*** argv,
    std::vector<int> &GlobalDim,
    std::vector<int> &GlobalPos
	){

	using namespace Grid;
	Grid_init(argc,argv);
	GridCartesian* UGrid = Grid::QCD::SpaceTimeGrid::makeFourDimGrid(
	GridDefaultLatt(), GridDefaultSimd(4, vComplex::Nsimd()),
	GridDefaultMpi());
	int NDIM=UGrid->_processors.size();

	GlobalDim.resize(NDIM);
	GlobalPos.resize(NDIM);
	for(int i=0;i<NDIM;i++){
		GlobalDim[i] = UGrid->_processors[i];
		GlobalPos[i] = UGrid->_processor_coor[i];
	}

	return NDIM;
	
}
#endif

int main (int argc, char** argv)
{
	
  std::vector<int> GlobalDim(4);
  std::vector<int> GlobalPos(4);

  int NDIM = init_QMP(&argc,&argv,GlobalDim,GlobalPos);
//  int NDIM = init_Grid(&argc,&argv,GlobalDim,GlobalPos);
  int loops;
  int sites = atoi (argv[1]); //global size
  int mem_size = atoi (argv[2]);
  int nblock = atoi (argv[3]);
  int verb = atoi (argv[4]);

    int GlobalIndex = BlockGeometry::CoorToIndex(GlobalPos,GlobalDim);

    BlockGeometry Src(NDIM);
    BlockGeometry Dest(NDIM);
	std::vector <int> TotalSites(NDIM); //number of sites

// Src is currently global, i.e. one for the whole
    for(int i=0;i<NDIM;i++){
		Src.BlockDim[i]=1;
		Dest.BlockDim[i]=1;
		TotalSites[i] = sites*GlobalDim[i];
	
	}	

	assert(nblock <=BlockGeometry::Total(GlobalDim));
	int temp_ind=0;
	while(nblock >1 ){
		if( Dest.BlockDim[temp_ind]  < GlobalDim[temp_ind] ){
		Dest.BlockDim[temp_ind] *=2;
		nblock = nblock/2;
		if(!GlobalIndex)printf("Dest.BlockDim[%d]=%d nblock=%d\n",temp_ind,Dest.BlockDim[temp_ind],nblock);
		}
		temp_ind = (temp_ind+1)%NDIM;
	}
//	std::cout << "Dest.BlockDim= "<<Dest.BlockDim;
//	std::cout <<std::endl;

	Src.SetLocal(GlobalPos,GlobalDim,TotalSites);
	Dest.SetLocal(GlobalPos,GlobalDim,TotalSites);


	if(!GlobalIndex) 
	printf ("GlobalDim = %d %d %d %d \n",
	GlobalDim[0], GlobalDim[1], GlobalDim[2], GlobalDim[3]);
	if(!GlobalIndex) 
	std::cout << "Dest.BlockDim= "<<Dest.BlockDim<<std::endl;

	if(0)
    PRINT ("GlobalPos = %d %d %d %d \n",
	GlobalPos[0], GlobalPos[1], GlobalPos[2], GlobalPos[3]);



    MPI_Comm mpi_comm; // used for global geometry
    MPI_Comm_split (MPI_COMM_WORLD,0,GlobalIndex,&mpi_comm);

    MPI_Comm mpi_split_comm;
    MPI_Comm_split (MPI_COMM_WORLD,Dest.BlockIndex(),Dest.NodeIndex(),&mpi_split_comm);


	if (Src.NodeTotal()!=(Dest.BlockTotal()*Dest.NodeTotal()))
	PRINT ("GlobalTotal=%d BlockTotal=%d LocalTotal=%d\n",
	1,Src.BlockTotal(),Src.NodeTotal());

//	Src.Print("Src");
//	Dest.Print("Dest");
	if(0)
	PRINT("Dest: Block (%d %d %d %d) %d Node (%d %d %d %d) %d\n",
	Dest.BlockPos[0],Dest.BlockPos[1],Dest.BlockPos[2],Dest.BlockPos[3],Dest.BlockIndex(),
	Dest.NodePos[0],Dest.NodePos[1],Dest.NodePos[2],Dest.NodePos[3],Dest.NodeIndex());
		

	double bytes = sizeof(DATA)*mem_size*Dest.DataVol();
	int VecTotal = Dest.BlockTotal();
	DATA * send_buf[VecTotal];
	for(int i=0;i<VecTotal;i++) 
		send_buf[i]= (DATA*)malloc(sizeof(DATA)*mem_size*Src.DataVol());
    DATA * recv_buf = (DATA*)malloc(sizeof(DATA)*mem_size*Dest.DataVol());
    DATA * recv2 = (DATA*)malloc(sizeof(DATA)*mem_size*Dest.DataVol());
	size_t offset1 = 10;
	while( offset1 <mem_size ) offset1 *=10;
	size_t offset2 = 10;
	while( offset2 <BlockGeometry::Total(TotalSites) ) offset2 *=10;
//	if(0)
	if (!GlobalIndex) std::cout << "offset1= "<<offset1<<" offset2= "<<offset2<<std::endl;

    for(size_t k=0;k<VecTotal;k++)
    for(size_t j=0;j<Src.DataVol();j++)
    for(size_t i=0;i<mem_size;i++){
		std::vector <int> SrcCoor(NDIM);
		BlockGeometry::IndexToCoor(j,SrcCoor,Src.DataDim);
		for(int dim=0;dim<NDIM;dim++){
			SrcCoor[dim] += Src.NodePos[dim]*Src.DataDim[dim];
		}
		size_t coor = BlockGeometry::CoorToIndex(SrcCoor,TotalSites);
	   *(send_buf[k] + i+mem_size*j )=i+offset1*(coor+offset2*k);
	}

    for(int j=0;j<Dest.DataVol();j++)
    for(int i=0;i<mem_size;i++){
	recv_buf[i+mem_size*j]=-1;
	recv2[i+mem_size*j]=-1;
	}

	std::vector<int> Index(VecTotal);
	std::vector<DATA *> sbuf(VecTotal);
	std::vector<DATA *> rbuf(1);
	rbuf[0] = recv_buf;
	for(int i=0;i<VecTotal;i++){
		Index[i] = i;
		sbuf[i] = send_buf[i];
	}
	Scramble<DATA> scr1(mem_size, &mpi_comm,verb);
	double t0 = dclock();
	scr1.run(Src,Index,sbuf,Dest,rbuf);
	double t1 = dclock();
	double bw = bytes/(t1-t0)/1000.; 
	if(!GlobalIndex) PRINT("scr1.run %g bytes / %g ms injection bw = %g MB/s per node \n",bytes,t1-t0,bw); t0=t1;


    for(size_t j=0;j<Dest.DataVol();j++)
    for(int i=0;i<mem_size;i++){
		std::vector<int> DestCoor(NDIM);
		BlockGeometry::IndexToCoor(j,DestCoor,Dest.DataDim);
		for(int k=0;k<NDIM;k++){
			DestCoor[k] +=Dest.NodePos[k]*Dest.DataDim[k];
		}
		size_t index = BlockGeometry::CoorToIndex(DestCoor,TotalSites);
		if( recv_buf[i+mem_size*j]!=(i+offset1*(index+offset2*Dest.BlockIndex())) ) 
    	PRINT ("recv_buf[%d][%d][%d] = %d (%d)\n",Dest.BlockIndex(),j,i,recv_buf[i+mem_size*j], (i+offset1*(index+offset2*Dest.BlockIndex())) );
    }
	t1 = dclock(); if(!GlobalIndex) PRINT("scr1.run check %g ms\n",t1-t0); t0=t1;


	Index.resize(1);
	sbuf.resize(1);
	rbuf.resize(VecTotal);
//    PRINT("Dest.BlockIndex()=%d\n",Dest.BlockIndex());
	Index[0] = Dest.BlockIndex();
	sbuf[0] = recv_buf;

	for(int i=0;i<VecTotal;i++){
		rbuf[i] = recv2+i*mem_size*Src.DataVol();
	}
	t1 = dclock(); if(!GlobalIndex) PRINT("scr1.run setup %g ms\n",t1-t0); t0=t1;
	scr1.run(Dest,Index,sbuf,Src,rbuf);
	t1 = dclock(); 
	bw = bytes/(t1-t0)/1000.; 
	if(!GlobalIndex) PRINT("scr1.run %g bytes / %g ms injection bw = %g MB/s per node \n",bytes,t1-t0,bw); t0=t1;
	
	if(1) 
    for(size_t k=0;k<Dest.BlockTotal();k++)
    for(size_t j=0;j<Src.DataVol();j++)
    for(size_t i=0;i<mem_size;i++)
	if( *(send_buf[k]+i+mem_size*j) != *(recv2+i+mem_size*(j+Src.DataVol()*k)) ) 
    {
    	PRINT ("send_buf[%ld][%ld][%ld] = %ld\n",k,j,i,*(send_buf[k]+i+mem_size*j));
    	PRINT ("recv2[%ld][%ld][%ld] = %ld\n",k,j,i,*(recv2+i+mem_size*(j+Src.DataVol()*k)));
    }


  
  std::cout <<scr1 <<"All passed!"<<std::endl;

  QMP_finalize_msg_passing ();

  return 0;
}

