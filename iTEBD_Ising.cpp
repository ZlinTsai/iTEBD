/**************************************************************************
Date: 2018/6/21

iTEBD algorithm for 1-dimension Ising model in transverse magnetic field.

- Zlin Tsai. s106022501@m106.nthu.edu.tw

***************************************************************************/

#include <iostream>
#include <sstream>
#include <fstream>
#include <uni10.hpp>
using namespace std;

uni10::UniTensor<double> SetBond(int v_in, int p_in, int v_out, int p_out, int D, int d)
{
    uni10::Bond virt_in(uni10::BD_IN, D);
    uni10::Bond phys_in(uni10::BD_IN, d);
    uni10::Bond virt_out(uni10::BD_OUT, D);
    uni10::Bond phys_out(uni10::BD_OUT, d);

    vector<uni10::Bond> bond;

    for (int i=0; i!=v_in; i++)
        bond.push_back(virt_in);

    for (int i=0; i!=p_in; i++)
        bond.push_back(phys_in);

    for (int i=0; i!=p_out; i++)
        bond.push_back(phys_out);

    for (int i=0; i!=v_out; i++)
        bond.push_back(virt_out);

    uni10::UniTensor<double> tensor(bond);

    return tensor;
}

uni10::UniTensor<double> OP(string name) 
{
	uni10::UniTensor<double> op = SetBond(0, 1, 0, 1, 1, 2); //cannot create a bond of dim = 0.

	if (name == "sz") 
    {
		op.SetElem({1.0, 0.0, \
					0.0, -1.0});
	}
	else if (name == "sx") 
    {
		op.SetElem({0.0, 1.0, \
					 1.0, 0.0});
	}
	else if (name == "id") 
    {
		op.SetElem({1.0, 0.0, \
					0.0, 1.0});
	}

	return op;
}

uni10::UniTensor<double> H_Ising(double h)
{
	uni10::UniTensor<double> sz = OP("sz");
	uni10::UniTensor<double> sx = OP("sx");
	uni10::UniTensor<double> id = OP("id");

    uni10::UniTensor<double> H = (-1.0) * uni10::Otimes(sz, sz) + (-0.5) * h * (uni10::Otimes(sx, id) + uni10::Otimes(id, sx));
    return H;
}

void SVD_MPS(uni10::UniTensor<double>& gammaA, uni10::UniTensor<double>& gammaB, \
             uni10::UniTensor<double>& lambdaA, uni10::UniTensor<double>& lambdaB, \
             uni10::UniTensor<double> theta, int D, int d, uni10::Network& svdL, uni10::Network& svdR) 
{
	vector<uni10::Matrix<double>> th_svd = uni10::Svd(theta.GetBlock()); //svd( )

	uni10::UniTensor<double> inv;
    inv = SetBond(1, 0, 1, 0, D, d);
    inv.PutBlock(uni10::Inverse(lambdaB.GetBlock()));
    svdL.PutTensor(0, inv);
    svdR.PutTensor(1, inv);

    //U
	uni10::Resize(th_svd[0], D*d, D, uni10::INPLACE);
	gammaA.PutBlock(th_svd[0]);
    svdL.PutTensor(1, gammaA);
    svdL.Launch(gammaA);

    //S
	uni10::Resize(th_svd[1], D, D, uni10::INPLACE);
	th_svd[1] *= (1.0/uni10::Norm(th_svd[1]));
    lambdaA.PutBlock(th_svd[1]);
    
    //V+
	uni10::Resize(th_svd[2], D, D*d, uni10::INPLACE); //careful have general from / D*d, d / D, D / D, D*d /
    uni10::Permute(gammaB, gammaB.label(), 1, uni10::INPLACE);
	gammaB.PutBlock(th_svd[2]);
    svdR.PutTensor(0, gammaB);
    svdR.Launch(gammaB);
}

double ExpValue(uni10::UniTensor<double>& gammaA, uni10::UniTensor<double>& gammaB, \
                uni10::UniTensor<double>& lambdaA, uni10::UniTensor<double>& lambdaB, \
                uni10::UniTensor<double> op, int idx, uni10::Network& state, uni10::Network& expVH, uni10::Network& expVS0, uni10::Network& expVS1) 
{
    uni10::UniTensor<double> expVal, ket;
    vector<uni10::UniTensor<double>> Vnorm, V;
    double norm; // <phi|phi>, uni10::Norm = sqrt(<phi|phi>) 

    if (idx == 0)
    { 
        Vnorm = {lambdaB, gammaA, lambdaA, gammaB, lambdaB};
        for (int i=0; i<Vnorm.size(); i++)
            state.PutTensor(i, Vnorm[i]);

        state.Launch(ket);
        norm = uni10::Norm(ket.GetBlock()) * uni10::Norm(ket.GetBlock());

        if (op.BondNum() == 4)
        {
            V = {lambdaB, gammaA, lambdaA, gammaB, lambdaB,\
                                    op,                    \
                 lambdaB, gammaA, lambdaA, gammaB, lambdaB};
            
            for (int i=0; i<V.size(); i++)
                expVH.PutTensor(i, V[i]);

            expVH.Launch(expVal);
            expVal *= (1.0/norm);
        }

        else if (op.BondNum() == 2)
        {
            V = {lambdaB, gammaA, lambdaA, gammaB, lambdaB,\
                            op,                            \
                 lambdaB, gammaA, lambdaA, gammaB, lambdaB};

            for (int i=0; i<V.size(); i++)
                expVS0.PutTensor(i, V[i]);
            
            expVS0.Launch(expVal);
            expVal *= (1.0/norm);
        }
	}
	else if (idx == 1)
    { 
        Vnorm = {lambdaA, gammaB, lambdaB, gammaA, lambdaA};
        for (int i=0; i<Vnorm.size(); i++)
            state.PutTensor(i, Vnorm[i]);

        state.Launch(ket);
        norm = uni10::Norm(ket.GetBlock()) * uni10::Norm(ket.GetBlock());

        if (op.BondNum() == 4)
        {
            V = {lambdaA, gammaB, lambdaB, gammaA, lambdaA,\
                                    op,                    \
                 lambdaA, gammaB, lambdaB, gammaA, lambdaA};
            
            for (int i=0; i<V.size(); i++)
                expVH.PutTensor(i, V[i]);
            
            expVH.Launch(expVal);
            expVal *= (1.0/norm);
        }

        else if (op.BondNum() == 2)
        {
            V = {lambdaB, gammaA, lambdaA, gammaB, lambdaB,\
                                             op,           \
                 lambdaB, gammaA, lambdaA, gammaB, lambdaB};
            
            for (int i=0; i<V.size(); i++)
                expVS1.PutTensor(i, V[i]);
            
            expVS1.Launch(expVal);
            expVal *= (1.0/norm);
        }
	}

	return expVal[0];
}

class iMPS
{
	private:

	public:
		iMPS(int, int);
		uni10::UniTensor<double> gammaA, gammaB, lambdaA, lambdaB, theta;
		//double ExpValue(uni10::UniTensor<double>, int, uni10::Network, uni10::Network, uni10::Network, uni10::Network);
		~iMPS() {};

};

iMPS::iMPS(int D, int d)
{
	gammaA = SetBond(1, 1, 1, 0, D, d);
    gammaB = SetBond(1, 1, 1, 0, D, d);
    lambdaA = SetBond(1, 0, 1, 0, D, d);
    lambdaB = SetBond(1, 0, 1, 0, D, d);

    //create random tensor
    gammaA.Randomize('U', -1, 1, 87);
    gammaB.Randomize('U', -1, 1, 87);
    uni10::Matrix<double> M( D, D, true); //set diagnoial random tensor
    M.Randomize('U', -1, 1, 114);
    lambdaA.PutBlock(M);
    lambdaB.PutBlock(M);

    uni10::Network state("network/State.net");
    uni10::Network svdL("network/SVD_L.net");
    uni10::Network svdR("network/SVD_R.net");

    vector<uni10::UniTensor<double>> V = {lambdaB, gammaA, lambdaA, gammaB, lambdaB};
    for (int i=0; i<V.size(); i++)
        state.PutTensor(i, V[i]);

    state.Launch(theta);
    SVD_MPS(gammaA, gammaB, lambdaA, lambdaB, theta, D, d, svdL, svdR); 
}



//iTEBD
uni10::UniTensor<double> Uevol(uni10::UniTensor<double> H, double dt)
{
	vector<uni10::Bond> bdH = H.bond();

	uni10::UniTensor<double> U(bdH);
	U.PutBlock(uni10::ExpH((-1.0)*dt, H.GetBlock()));

	return U;
}

void updateMPS(uni10::UniTensor<double>& gammaA, uni10::UniTensor<double>& gammaB, uni10::UniTensor<double>& lambdaA, uni10::UniTensor<double>& lambdaB,
			   uni10::UniTensor<double> U, int D, uni10::Network& stateU, uni10::Network& svdL, uni10::Network& svdR)
{
	int d = U.bond()[0].dim();

    vector<uni10::UniTensor<double>> V = {lambdaB, gammaA, lambdaA, gammaB, lambdaB,\
                                                              U                    };
    for (int i=0; i<V.size(); i++)
        stateU.PutTensor(i, V[i]);
    
	uni10::UniTensor<double> theta;
    stateU.Launch(theta);
    SVD_MPS(gammaA, gammaB, lambdaA, lambdaB, theta, D, d, svdL, svdR);
}

void trotterSuzuki(uni10::UniTensor<double>& gammaA, uni10::UniTensor<double>& gammaB, uni10::UniTensor<double>& lambdaA, uni10::UniTensor<double>& lambdaB,
				   uni10::UniTensor<double> U, int D, uni10::Network& stateU, uni10::Network& svdL, uni10::Network& svdR) 
{
	updateMPS(gammaA, gammaB, lambdaA, lambdaB, U, D, stateU, svdL, svdR);
	updateMPS(gammaB, gammaA, lambdaB, lambdaA, U, D, stateU, svdL, svdR);
}

class iTEBD_1D
{
	private:

	public:
		iTEBD_1D(iMPS&, uni10::UniTensor<double>, int, double);
};

iTEBD_1D::iTEBD_1D(iMPS& param, uni10::UniTensor<double> H, int D, double dt)
{
    uni10::Network state("network/State.net");
    uni10::Network stateU("network/Theta.net");
    uni10::Network svdL("network/SVD_L.net");
    uni10::Network svdR("network/SVD_R.net");
    uni10::Network expVH("network/ExpValueH.net");
    uni10::Network expVS0("network/ExpValueS0.net");
    uni10::Network expVS1("network/ExpValueS1.net");
    
	uni10::UniTensor<double> U;
	double E_new = 1.0, E_old = 0.0, diff = 1.0, times = 0.0; 
	int x = -8, steps = 0;
	for (int i=1; i<5; i++)
    {
		dt = dt/2.0;
		U = Uevol(H, dt);
		while (diff >= pow(10, x))
        {
			steps++;
			//times += dt;
			trotterSuzuki(param.gammaA, param.gammaB, param.lambdaA, param.lambdaB, U, D, stateU, svdL, svdR);
            if (steps%100 == 0 || steps%101 == 0 || steps >= 10000)
            {
                E_old = E_new;
                E_new = 0.5 * (ExpValue(param.gammaA, param.gammaB, param.lambdaA, param.lambdaB, H, 0, state, expVH, expVS0, expVS1) + \
                               ExpValue(param.gammaA, param.gammaB, param.lambdaA, param.lambdaB, H, 1, state, expVH, expVS0, expVS1));
                diff = abs(E_old - E_new);
            }
		}
		x -= 1;

		if (diff <= pow(10, -10)) //conv condition -10 is good for D = 5~50, D = 100 maybe need to -12
			break;
	}
}

int main()
{
    double h = 1.0, dt = 0.2, f;
	uni10::UniTensor<double> H = H_Ising(h);
	int d = H.bond()[0].dim();
	int D = 20;

    iMPS State(D, d);
    vector<double> hs, En, Sz0, Sz1, Sx0, Sx1;
    
	for(f=0; f<=15; f+=1) // change h field  for(f=1000;f<=1000;f+=2)
	{
		h = 0.1*f; 	
		H = H_Ising(h);
    
        iTEBD_1D iTEBD(State, H, D, dt);

        uni10::Network state("network/State.net");
        uni10::Network expVH("network/ExpValueH.net");
        uni10::Network expVS0("network/ExpValueS0.net");
        uni10::Network expVS1("network/ExpValueS1.net");
    
        double en = 0.5 * (ExpValue(State.gammaA, State.gammaB, State.lambdaA, State.lambdaB, H, 0, state, expVH, expVS0, expVS1) + \
                           ExpValue(State.gammaA, State.gammaB, State.lambdaA, State.lambdaB, H, 1, state, expVH, expVS0, expVS1));
        double sz0 = ExpValue(State.gammaA, State.gammaB, State.lambdaA, State.lambdaB, OP("sz"), 0, state, expVH, expVS0, expVS1);
        double sz1 = ExpValue(State.gammaA, State.gammaB, State.lambdaA, State.lambdaB, OP("sz"), 1, state, expVH, expVS0, expVS1);
        double sx0 = ExpValue(State.gammaA, State.gammaB, State.lambdaA, State.lambdaB, OP("sx"), 0, state, expVH, expVS0, expVS1);
        double sx1 = ExpValue(State.gammaA, State.gammaB, State.lambdaA, State.lambdaB, OP("sx"), 1, state, expVH, expVS0, expVS1);

        hs.push_back(h); En.push_back(en); Sz0.push_back(sz0); Sz1.push_back(sz1); Sx0.push_back(sx0); Sx1.push_back(sx1); 
        cout << setprecision(8) << h << " " << en << " " << sz0 << " " << sz1 << " " << sx0 << " " << sx1 << endl;
    }

    // save data
    string file;
    file = "database/bondDim_" + to_string(D) + ".csv";
    cout << "\nThe data already save at " << file << endl;
    ofstream fout(file);
    if (!fout)
        cout << "Fail to save (maybe need mkdir data/)" << endl;

    fout << "h,En,Sz_0,Sz_1,Sx_0,Sx_1" << endl;
    fout.flush();
    for (int i=0; i<En.size(); i++)
    {
        fout << hs[i] << "," << En[i] << "," << Sz0[i] << "," << Sz1[i] << "," << Sx0[i] << "," << Sx1[i] << endl;
    }
    fout.flush();
    fout.close();

    return 0;
}
