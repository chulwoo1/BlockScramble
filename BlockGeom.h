#ifndef BLOCK_GEOM_H
#define BLOCK_GEOM_H

template<typename T>
std::ostream& operator<<(std::ostream&s, std::vector<T> t) {
		s<<"[";
	for(auto i: t) { 
		s << " " <<t[i] ;
	}
		s<<"]"<<std::endl;
}

class 	CartesianGeometry {
public:
static	size_t Total(std::vector <int > & Dim){
	unsigned int ndim = Dim.size();
	size_t tmp=1;
	for(int i = 0;i<ndim;i++){
		tmp = tmp*Dim[i];
	}
	return tmp;
	}

static	size_t CoorToIndex(std::vector <int> & Pos,std::vector <int > & Dim){
//	std::cout <<"Pos "<<Pos;
//	std::cout <<"Dim "<<Dim;
	unsigned int ndim = Pos.size();
	assert(ndim==Dim.size());
	int tmp = Pos[ndim-1];
	for(int i = ndim-2;i>=0;i--){
		tmp = tmp*Dim[i] + Pos[i];
	}
	return tmp;
	}

static	int IndexToCoor( size_t index, 
				std::vector <int> & Pos,
				std::vector <int > & Dim){
	unsigned int ndim = Pos.size();
	assert(ndim==Dim.size());
	int tmp = index;
	for(int i = 0;i<ndim;i++){
		Pos[i] = tmp % Dim[i];
		tmp = tmp / Dim[i];
	}
	return tmp;
	}

};

class BlockGeometry : public virtual CartesianGeometry {
	public:
	std::vector<int>  NodeDim;
	std::vector<int>  NodePos;
	std::vector<int>  BlockDim;
	std::vector<int>  BlockPos;
	std::vector<int>  DataDim;
	int ndim;
	BlockGeometry() {}

	BlockGeometry(int dim):ndim(dim){ resize(dim);}
	int Dim(){ return ndim;}
   	int resize(int dim){
		NodeDim.resize(dim);
		NodePos.resize(dim);
		BlockDim.resize(dim);
		BlockPos.resize(dim);
		DataDim.resize(dim);
		ndim=dim;
//		std::cout <<"BlockGeometry(): ndim= "<<ndim<<std::endl;
		return ndim;
	}
	void Check(std::vector <int> &Dim){ assert(Dim.size()==ndim);};

	void SetLocal(std::vector <int> & GlobalPos,std::vector <int > & GlobalDim,
			std::vector<int> & TotalSites){
		Check(GlobalPos);Check(GlobalDim);
		Check(TotalSites);
		for(int i=0;i<ndim;i++){
//			std::cout <<" GlobalDim "<<i<<" "<<GlobalDim[i];
//			std::cout <<" GlobalPos "<<i<<" "<<GlobalPos[i];
//			std::cout <<" BlockDim "<<i<<" "<<BlockDim[i];
//			std::cout <<std::endl;
			assert(GlobalDim[i]%BlockDim[i]==0);
			NodeDim[i] = GlobalDim[i]/BlockDim[i];
			NodePos[i] = GlobalPos[i]%NodeDim[i];
			BlockPos[i] = GlobalPos[i]/NodeDim[i];
			DataDim[i] = TotalSites[i]/NodeDim[i];
		}
	}

	size_t BlockTotal(){ return Total(BlockDim); }
	size_t NodeTotal(){ return Total(NodeDim); }
	size_t DataVol(){ return Total(DataDim); }
	size_t BlockIndex(){ return CoorToIndex(BlockPos,BlockDim);}
	size_t NodeIndex(){ return CoorToIndex(NodePos,NodeDim);}

	void Print(std::string name){

	std::cout << name <<" NodeDim:";
	for(auto i : NodeDim) { std::cout << ' '<< i; }
	std::cout<<std::endl;
	std::cout << name <<" DataDim:";
	for(auto i : DataDim) { std::cout << ' '<< i; }
	std::cout<<std::endl;
	}

};


#endif
