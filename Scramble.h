#ifndef SCRAMBLE_H
#define SCRAMBLE_H

#include <mpi.h>
//#include <qmp.h>
#include "BlockGeom.h"

template < typename DATA > 
class Scramble : public virtual CartesianGeometry{
public:

  size_t mem_size;
  MPI_Comm *mpi_comm;
  int verb;
  int rank, size;

  Scramble (int _mem_size, MPI_Comm * _mpi_comm, int _verb=0)
	:mem_size (_mem_size), mpi_comm (_mpi_comm),verb(_verb)
  {
	MPI_Comm_rank(*mpi_comm,&rank);
	MPI_Comm_size(*mpi_comm,&size);
  }

  void SetVerbosity(int _verb) {verb=_verb;}

#if 1
  friend std::ostream& operator<<(std::ostream&s, Scramble &scr ) {
	s<< "Scramble: " << scr.rank <<" of "<<scr.size<<" : ";
        return s;
  }
#endif

  int GCF(int a, int b){

	int small,large;
	if (a>b){ large=a;small=b;}
	else { large=b;small=a;}
	if (verb>6)
	std::cout << "GCF: "<<large<<" "<<small<<std::endl;
	
	int tmp = (large%small);
	while( tmp > 0 ){
		std::cout << "GCF: "<<large<<" "<<small<<std::endl;
		large=small;
		small=tmp;
		tmp = (large%small);
	}

	if (small<1) small=1;

	return small;
	

  }
 

  void run (BlockGeometry & Src,
	    std::vector < int >&SrcIndex,
	    std::vector < DATA * >send_buf, BlockGeometry & Dest,
	    std::vector < DATA * >recv_buf)
  {

    assert (SrcIndex.size () == send_buf.size ());
	int bsize = GCF(Src.DataDim[0],Dest.DataDim[0]);
//	bsize=1;
	if (verb>5 && !rank )
	std::cout << "bsize: "<<bsize<<std::endl;
    int NDIM = Src.Dim ();
    int recv_max = recv_buf.size ();
for (int recv_i = 0; recv_i < recv_max; recv_i++) {
      MPI_Win  recv_win;
	MPI_Win_create (recv_buf[recv_i], sizeof (DATA) * mem_size * Dest.DataVol (),
			sizeof (DATA)*mem_size, MPI_INFO_NULL, *mpi_comm,
			&recv_win);
	MPI_Win_fence (MPI_MODE_NOPRECEDE, recv_win);

      for (int k = 0; k < SrcIndex.size (); k++) {
	  std::vector < int >DestBlockCoor (NDIM);
	  IndexToCoor (SrcIndex[k], DestBlockCoor, Dest.BlockDim);
	  int DestWrap = SrcIndex[k]/Dest.BlockTotal();

	if (DestWrap != recv_i) continue;
	if(verb>5)
	printf("rank %d: SrcIndex[%d]=%d DestWrap=%d recv_i=%d\n",rank,
	k,SrcIndex[k],DestWrap,recv_i);

	if(!rank)
	if (verb>5)
      std::cout<<"SrcIndex "<<SrcIndex[k]<<" Src.BlockTotal "<<Src.BlockTotal()<< " Src.BlockIndex "<<Src.BlockIndex()<<std::endl;
	// check to see if the SrcIndex is eligible for the block
	assert ((SrcIndex[k] % Src.BlockTotal ()) == Src.BlockIndex ());
	std::vector < int >GlobalDim (NDIM);
	for (int i = 0; i < NDIM; i++) {
	  GlobalDim[i] = Dest.BlockDim[i] * Dest.NodeDim[i];
	  assert (GlobalDim[i] == (Src.BlockDim[i] * Src.NodeDim[i]));
	}


#pragma omp parallel for 
	for (size_t j = 0; j < Src.DataVol (); j += bsize) {
	if(verb>5)
	std::cout <<j<<" : thread "<<omp_get_thread_num()<<" of "<<omp_get_num_threads()<<std::endl;
	  std::vector < int >SrcCoor (NDIM);	//global site coordinate
	  IndexToCoor (j, SrcCoor, Src.DataDim);
	  for (int i = 0; i < NDIM; i++) {
	    SrcCoor[i] += Src.NodePos[i] * Src.DataDim[i];
	  }


	  //In case destintation receives multiple vectors, currrently assumed to be contiguous

	  std::vector < int >DestNodeCoor (NDIM);
	  std::vector < int >DestSiteCoor (NDIM);
	  for (int i = 0; i < NDIM; i++) {
	    DestSiteCoor[i] = SrcCoor[i] % Dest.DataDim[i];
	    DestNodeCoor[i] = SrcCoor[i] / Dest.DataDim[i];
	    DestNodeCoor[i] += DestBlockCoor[i] * Dest.NodeDim[i];
	  }


	  size_t target = CoorToIndex (DestNodeCoor, GlobalDim);
	  size_t offset = CoorToIndex (DestSiteCoor, Dest.DataDim);
	if(verb>3 &&  !rank)
	    if (target == 0 && offset == 0) {
	      printf
		("target 0 offset 0 DestWrap %d Index %d DestNodeCoor %d %d %d %d DesSiteCoor %d %d %d %d\n",
		 DestWrap, SrcIndex[k], DestNodeCoor[0], DestNodeCoor[1],
		 DestNodeCoor[2], DestNodeCoor[3], DestSiteCoor[0],
		 DestSiteCoor[1], DestSiteCoor[2], DestSiteCoor[3]);
	    }
	  assert (DestWrap < recv_max);

#pragma omp critical
{
	  MPI_Put (send_buf[k] + mem_size * j, bsize*mem_size * sizeof (DATA),
		   MPI_BYTE, target, offset,
		   bsize*mem_size * sizeof (DATA), MPI_BYTE, recv_win);
	if(verb>6)
	std::cout << *this << "MPI_Put: send_buf["<< k << "]+" <<mem_size * j<<" "<<bsize*mem_size * sizeof (DATA) << " : "<<target <<" "<< offset <<" "<<bsize*mem_size * sizeof (DATA) <<std::endl;
}
#if 0
	  printf ("rank %d: MPI_Put: send_buf[%d]+%d %d : %d %d %d\n",
		       rank, k, mem_size * j, bsize*mem_size * sizeof (DATA),
		       target, offset , bsize*mem_size * sizeof (DATA));
#endif
	}
      }
	MPI_Win_fence ((MPI_MODE_NOSTORE | MPI_MODE_NOSUCCEED), recv_win);
	MPI_Win_free (&recv_win);
}
  }

};

#endif
