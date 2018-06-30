/*

phi4.cpp
Definitions declared in phi4.h
for a massive scalar Phi^4 simulation using the Truncated Hilbert Space (THS) method

Author: Marton Lajer, Physics MSC student
Supervisor: Zoltan Bajnok, Research Professor

Institute for Theoretical Physics
Roland Eotvos University, Budapest (2014-2015)

For more information, please read the heading of phi4.h!

*/


#include "phi4.h"

#include <iostream>

double sqr(double d)
{
	return d*d;
}

bool compare(elem e1, elem e2)				//  elements in a HilbertState will be sorted firstly by |k|, then negative momenta first: eg. 0,-1,+1,-2,+2,...
{
	if (abs(e1.k) != abs(e2.k))
		return (abs(e1.k) < abs(e2.k));
	return e1.k < e2.k;
}

bool HilbertSpace::HScomp(HilbertState H1, HilbertState H2)
{
        return (Energy0(H1) < Energy0(H2));
}

//----------------------[Methods for HilbertState class]----------------------------//


HilbertState::HilbertState()
{
	particle_num=0;
}

HilbertState::HilbertState(std::vector<elem> state_in)
{
	particles = state_in;
	particle_num=0;
	std::sort(particles.begin(),particles.end(),compare);
	for (int i=0; i < particles.size(); i++)
		particle_num += particles[i].n;
}

double HilbertSpace::Energy0(const HilbertState& HS)
{
	if (HS.StateNum() == 0) return 0.0;
	double E0=0.0;
	if (periodicBC==1)
	{
        for (unsigned int i=0; i<HS.StateNum(); i++)
        {
            elem e = HS.getState(i);
            E0+=(double)e.n*sqrt(sqr(M)+sqr((double)e.k*2.0*pi/L));
        }
     } else if (periodicBC == -1)
     {
        for (unsigned int i=0; i<HS.StateNum(); i++)
        {
            elem e = HS.getState(i);
            if (e.k > 0)
            {
                E0+=(double)e.n*sqrt(sqr(M)+sqr(((double)e.k*2.0-1.0)*pi/L));
            }
            else
            {
                E0+=(double)e.n*sqrt(sqr(M)+sqr(((double)e.k*2.0+1.0)*pi/L));
            }
        }
     }
	return E0;
}

int HilbertState::ParticleNumber() const
{
	return particle_num;
}

int HilbertState::max_k()
{
	if (particles.size() == 0) return 0;
	return particles[particles.size()-1].k;
}

int HilbertState::Create(int k, int n)
{
		for (unsigned int i=0; i < particles.size(); i++)
		{
			if (particles[i].k == k)
			{
				if ((-n) > particles[i].n) { particle_num -= particles[i].n; particles[i].n=0; return -1; }
				particles[i].n += n;
				particle_num += n;
				return particles[i].n;
			}
		}
		if (n >= 0)
		{
			elem e = elem(k,n);
			particles.push_back(e);
			particle_num += n;
			std::sort(particles.begin(),particles.end(),compare);
		}
		return n;
}


//----------------------------[Methods for HilbertSpace class]----------------------------//


HilbertSpace::HilbertSpace()
{
}

HilbertSpace::HilbertSpace(std::vector<HilbertState> space_in)
{
	states = space_in;
}

void HilbertSpace::MomentumZero(HilbertState state_in, int PBC_in, short momSec)
{
	int iterator = 0;
	HilbertState neg;
	for (int j = 0; j < state_in.StateNum(); j++)
	{
		elem e = state_in.getState(j);
		neg.Create(-e.k,0);
	}
	bool flag = true;
	do
	{
		if (flag)
		{
		int SUM =0;
		for (int j = 0; j < state_in.StateNum(); j++)
		{

			elem e1 = state_in.getState(j);
			elem e2 = neg.getState(j);
            if (PBC_in == 1)
            {
                SUM += e1.k*(e1.n-2*e2.n);
            } else
            if (PBC_in == -1)
            {
                if (e1.k > 0)
                {
                    SUM += (2*e1.k-1)*(e1.n-2*e2.n);
                } else
                {
                    SUM += (2*e1.k+1)*(e1.n-2*e2.n);
                }
            }
		}

		if (SUM == momSec)   // SUM == 0
		{
			HilbertState hs;
			for (int j = 0; j < state_in.StateNum(); j++)
			{
				elem e1 = state_in.getState(j);
				elem e2 = neg.getState(j);
				if (e2.n > 0)
					hs.Create(-e1.k,e2.n);
				if ((e1.n-e2.n)>0)
					hs.Create(e1.k,e1.n-e2.n);
			}
			if (states.size() == 0)
			{
				states.push_back(hs);
			} else
			if (hs != states[states.size()-1])
            {
			    states.push_back(hs);
			}
		}
		}


		elem e= state_in.getState(iterator);
		elem f= neg.getState(iterator);

		if (f.n < e.n)
		{
			neg.Create(-e.k,1);
			iterator=0;
			flag = true;
		} else
		if (f.n == e.n)
		{
			iterator++;
			neg.Create(-e.k,-e.n);
			flag = false;
		}

	} while (iterator < state_in.StateNum());
}

bool HilbertSpace::FillStates(HilbertState start, short nonMiniOnly, double Energy0Max, short parity, short momSec)
{
	int i=start.max_k();
	if ((i==0)&&((periodicBC==-1)||nonMiniOnly))
        i = 1;  //nonzero momenta only!
	do
	{
		HilbertState v = start;
		v.Create(i,1);
		if (Energy0(v) >= Energy0Max)
			return false;
		FillStates(v,nonMiniOnly,Energy0Max,parity, momSec);
		if ((parity==2)||((v.ParticleNumber()%2)==parity))
            MomentumZero(v,periodicBC, momSec);
		i++;
	} while (1);
}

int HilbertSpace::NumStates() const
{
	return states.size();
}


void HilbertSpace::DimTrunc(unsigned int prefDim)
{
    std::sort(states.begin(),states.end(),HScomp);
    double Emax = Energy0(prefDim-1);
    unsigned int i =prefDim-1;
    while (Energy0(i)-Emax < 1e-12)
    {
        i++;
        if (i==states.size()) break;
    }
    if (i < states.size())
        states.resize(i);
}

void HilbertSpace::DimTrunc_nosort(unsigned int prefDim)
{
    double Emax = Energy0(prefDim-1);
    unsigned int i =prefDim-1;
    while (Energy0(i)-Emax < 1e-12)
    {
        i++;
        if (i==states.size()) break;
    }
    if (i < states.size())
        states.resize(i);
}

int HilbertSpace::FindState(HilbertState& hs_in, int max_ind)
{
    double e0=Energy0(hs_in);
    int ind0 = 0;
    int ind1 = max_ind;
    int ind = 0;
    do
    {
        ind = ind0+(ind1-ind0)/2;

        if (Energy0(ind)>(e0+1e-12))
            ind1 = ind;
        else
        if (Energy0(ind)<(e0-1e-12))
            ind0 = ind;

        if (ind1-ind0 <= 1)
        {
            if (states[ind0]==hs_in)
                return ind0;
            if (states[ind1]==hs_in)
                return ind1;
            return -1;
        }
    }
    while (std::abs(Energy0(ind)-e0)>1e-12);
    ind1=ind;
    ind0=ind1;
    while ((Energy0(ind0) > e0-1e-12)&&(ind0>0))
        ind0--;
    while ((Energy0(ind1) < e0+1e-12)&&(ind1<max_ind))
        ind1++;

    for (ind = ind0; ind <= ind1; ind++)
    {
        if (states[ind]==hs_in)
            return ind;
    }


    return -1;
}

bool indsort(std::vector<int> a, std::vector<int> b)
{
    return a[a.size()-1]<b[b.size()-1];
}

bool intsort(int a, int b)
{
    return a<b;
}

//----------------------------------------------------------------------------------------//

double basis(HilbertSpace &H_in, HilbertState vacuum_in, short useMiniSpace, double E0, int preferredDimension, short parity, short momSec)
{
    short par = parity;
    if (useMiniSpace)
        par = 2;
    HilbertSpace H;
    if (preferredDimension == 0)
    {
        H_in = H;
         return 0;
    }

    H.FillStates(vacuum_in, useMiniSpace, E0, par, momSec);
    int n=H.NumStates();

    if (H.NumStates()<preferredDimension)
    {
        H.clear();
        H.FillStates(vacuum_in,useMiniSpace, E0*3/4, par, momSec);

        int n1=H.NumStates();
        if (n1==0)
        {
            std::cerr << "Please give a bigger energy guess!\n";
            return -1;
        }

        double Energy = E0*(1.0+0.25*log((double)preferredDimension*1.1/n)/log((double)n/n1));
        H.clear();
        H.FillStates(vacuum_in,useMiniSpace, Energy, par, momSec);
        if (H.NumStates() < preferredDimension)
            basis(H_in,vacuum_in,useMiniSpace, Energy,preferredDimension, par, momSec);
    }
    if (H.NumStates() >= preferredDimension)
    {
        if ((momSec==0)&&((par == 0)||(par==2)))
            H.AddState(vacuum_in);
        H.DimTrunc(preferredDimension);
        H_in = H;
    }

    return H_in.Energy0(H_in.NumStates()-1);
}
