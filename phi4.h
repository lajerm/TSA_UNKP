/*

phi4.h
Declarations and basic definitions for classes HilbertState and HilbertSpace

Author: Marton Lajer, Physics MSC student
Supervisor: Zoltan Bajnok, Research Professor

Institute for Theoretical Physics
Roland Eotvos University, Budapest (2014)

*/


#ifndef __PHI_4_H__
#define __PHI_4_H__

#include <vector>
#include <cmath>
#include <algorithm>

const double pi = 3.14159265358979;

double sqr(double d);

struct elem
/*
	The basic building block of a HilbertState. Can store a pair of integers:
	k: numerical momentum of a particle (e.g. in the periodic sector, it's values are  k*L/(2*pi) ->  k=...,-2,-1,0,1,2,... )
	n: number of particles with k numerical momentum, -> n=0,1,2,...
*/
{
	int k;
	int n;

    // constructors
	elem() {k=0; n=0;}
	elem(int k_in, int n_in) {k=k_in; n=n_in;}
};

bool compare(elem e1, elem e2);

class HilbertState
/*
	This class realizes an eigenstate of the free Hamiltonian (an element of the Fock space).
	It is structured around a vector<elem>.
*/
{
	std::vector<elem> particles;
	int particle_num;

public:

//constructors:
	HilbertState();
	// creates an "empty" state (the vacuum):

	HilbertState(std::vector<elem>);
	// creates a state with several particles (see -> elem)

// destructor:
	~HilbertState() {particles.clear();}

// member functions:
	int ParticleNumber() const;
	//Gives the particle number of the state

	int max_k();
	// Gives the numerical momentum  (see -> elem) of the particle(s)  in the state that have the largest absolute value of k

	int Create(int k, int n);  //n < 0-ra annihilál!
	// Creates (or annihilates!) n particles with k numerical momentum. Returns the new number of k-momentum particles in the state.
    // If (-n) > the number of k-momentum particles -> returns (-1) and annihilates all the k-momentum particles in the state.
    // (for an explanation, see the calculation of matrix elements in phi_n_interaction::n_MxEl(...) in phi_n.cpp)

	elem getState(int num) const { return particles[num]; }
	// Gives the state stored in the 'num'th. element of the underlying private vector<elem> 'particles'.

	void setState(int num, elem e_in) { particle_num += e_in.n-particles[num].n; particles[num] = e_in; return; }
	// Changes the state stored in the 'num'th. element of the vector 'particles' to e_in.

	int StateNum() const { return particles.size(); }
	// Returns the size of the vector 'particles'.

// overloaded operators:
	bool operator != (const HilbertState& other) const
	// Test whether two HilbertStates are different
	{
		if (particles.size() != other.particles.size()) return true;
		for (unsigned int j=0; j < particles.size(); j++)
		{
			if ((particles[j].k !=other.particles[j].k)||(particles[j].n !=other.particles[j].n))
				return true;
		}
		return false;
	}

    // Test if two HilbertStates are the same
	bool operator == (const HilbertState& other) const
	{
	    return (!(*this!=other));
	}

	void killZeros()
	{
	    int j=0;
	    for (int i=0; i < particles.size(); i++)
        {

            while ((particles[i+j].n==0)&&(i+j < particles.size()))
            {
                j++;
            }
            if (i+j < particles.size())
            {
                particles[i]=particles[i+j];
            }
            else
                break;
        }
        particles.resize(particles.size()-j);
	}

};

class HilbertSpace
/*
    This class realizes the truncated Hilbert space in which we are to diagonalize the Hamiltonian.
	Its data structure is a vector<HilbertState>.
*/
{
	std::vector<HilbertState> states;
	static double L;
	static double M;
	static int periodicBC;

public:
// Constructors:
	HilbertSpace();
	// constructs an empty variable

	HilbertSpace(std::vector<HilbertState>);
	// constructs the Hilbert space with predefined base elements.

// Destructor:
	~HilbertSpace() {states.clear(); }

// Member functions:
	bool FillStates(HilbertState start, short useMiniSpace, double Energy0Max, short parity, short momSec);
	// Recursively generates all the Fock vectors below a given Energy0Max,
    // in the given particle parity and momentum sector. useMiniSpace is set to 1
    // if the minimal space method is used, otherwise 0.

	void MomentumZero(HilbertState state_in, int PBC_in, short momSec);
	// auxiliary function that assists FillStates(...). FillStates calls MomentumZero
    // with input states where all particles have positive momentum. From such a state,
    // MomentumZero calculates all states with the same |k|'s (and thus with the same energy),
    // but with p=0. These states are then written to the vector<HilbertState>.

    void nonzero(int index0, std::vector<int> &indices);        // early development version for 2nd and 4th order interaction only

    int FindState(HilbertState& hs_in, int max_ind);

	void AddState(HilbertState state_in) {states.push_back(state_in); }
	// Explicitly adds an element to the end of the vector<HilbertState>.

	void clear() {states.clear(); }
	// Clears the vector<HilbertState>.

	int NumStates() const;
	//  Gives the dimension of the truncated Hilbert space.

	HilbertState GetState(int n) const {return states[n];}
    // Gives the 'n'th base element stored in the vector<HilbertState>.

	static double Energy0(const HilbertState& HS);
	// Calculate the free Hamiltonian eigenvalue corresponding to state HS

	double Energy0(int i) const {return Energy0(states[i]);}
	// Same as the previous, i is the index of the HilbertState in the vector 'states'

	void DimTrunc(unsigned int prefDim);
	// Truncate the basis at a given dimension number prefDim

	void DimTrunc_nosort(unsigned int prefDim);
    // same as the previous but the states are not sorted

	static void setParams(double L_in, double M_in, int PBC_in) {L=L_in; M=M_in; periodicBC = PBC_in;}
    // set the static parameters L, M and PBC

	static double getL() {return L;}
	static double getM() {return M;}
	static int getBC() { return periodicBC;}

	static bool HScomp(HilbertState H1, HilbertState H2);
	// compare two eigenstates with respect to their free Hamiltonian eigenvalues. Returns true if
	// Energy0(H1) < Energy0(H2)
};

double basis(HilbertSpace &H_in, HilbertState vacuum_in, short useMiniSpace, double E0, int preferredDimension, short parity, short momSec);
// Creates a basis with given preferred dimension in a given parity and overall momentum sector, and stores it in
// the HilbertSpace variable H_in. Returns with the free Ham. energy of the highest energy state.
// The method needs an initial energy guess E0 for the energy cutoff


#endif
