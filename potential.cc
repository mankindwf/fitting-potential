#include <fstream>       
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <cassert>
#include <vector>
#include <cmath>
#include <mpi.h>
#include <infer/state.h>
#include <infer/objfunction.h>
#include <infer/pso.h>
#define LBOX 3.61*3

using namespace std;
using namespace pso;

typedef std::vector<double> config;

class Potential{
public:
	virtual double operator()(const config &conf, const State &param)=0;
	int N;
};


double PairPotential(double r, double epsilon, double sigma){
	return 4.0*epsilon*(pow(sigma/r,12) - pow(sigma/r,6) );}



double Energia(const config &conf, const State &param){
		double LHALF = 0.5*LBOX;
		int natoms = conf.size()/3;
		double U = 0.0e0;
		for(int i=0;i<natoms-1;++i){for(int j=i+1;j<natoms;++j){
    		double dr[3], r2=0.0e0;
    		for (int q=0;q<3;++q) dr[q] = conf[j*3+q]- conf[i*3+q];
    		for (int q=0;q<3;++q){
     		double * dd = dr+q;
     		if (*dd >= LHALF) (*dd) -= LBOX;
     		else if (*dd < -LHALF) (*dd) += LBOX;
     		r2 += (dr[q]*dr[q]);
    		}
    		U += PairPotential(sqrt(r2), param[0], param[1]);
    		}}
 	return U;}


class LJ: public Potential{
public:
	LJ() { N = 2; }
	double operator()(const config &conf, const State &param){return Energia(conf, param);}
};



void LeerPos(char *data_pos, vector <config> &all_atoms){
		ifstream dataxyz;
		dataxyz.open(data_pos);
		int i=0;
		while( 1 ){
			int num_atoms;
			dataxyz >> num_atoms;
			if( dataxyz.eof() ){break;}
			//cout << num_atoms << endl;
			//cout << "TEST" << endl;
			config atoms(3*num_atoms);
		for(int k=0; k<num_atoms; k++){dataxyz >> atoms[3*k] >> atoms[3*k+1] >> atoms[3*k+2];}
			all_atoms.push_back(atoms);
			i++;}}


class Error: public ObjectiveFunction{
public:
	Error(char *data_pos, char *data_ene, Potential &Pot): ObjectiveFunction(Pot.N), potential(Pot)
               {

		LeerPos(data_pos, all_atoms);
		ifstream dataene;
		dataene.open(data_ene);
		double each_energy;
		while(!dataene.eof()){  dataene >> each_energy;
					energies.push_back(each_energy);}
		}

	double operator()(const State &param)const{	
		double total=0.0;		
		int num_config = all_atoms.size();
		for(int i=0; i<num_config; i++){
				total += pow(potential(all_atoms[i], param) -energies[i],2);}
		return total;
	}

private:
	vector <config> all_atoms;
	config energies;
        Potential & potential;
	
};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//					 	INICIO	 						////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv){
MPI_Init(&argc, &argv);


vector <config> all_atoms;
LeerPos("atoms.xyz", all_atoms);
int num_config = all_atoms.size();

LJ Pot;
State param({1,3.4});
for(int i=0; i<num_config; i++){cout << Pot(all_atoms[i], param) << endl;}

Error MinLJ("atoms.xyz","energias.dat",Pot);
State s({1.5,3.5});
cout << MinLJ(s) << endl;

Minimizer min(2, 100);
min.Minimize(MinLJ, s, 1e-6);

MPI_Finalize();
return 0;}
