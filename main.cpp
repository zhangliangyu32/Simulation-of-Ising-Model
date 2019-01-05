#include <iostream>
#include <cmath>
#include <ctime>
#include <fstream>
#include <cstdlib>
#include <iostream>
using namespace std;
class State
{
private:
    int * configuration;
    int M;//共有M * M个节点
    double beta;
    double J;
    double h;
public:
    State(int M0, double beta0, double J0, double h0)
    {
        M = M0;
        beta = beta0;
        J = J0;
        h = h0;
        configuration = new int[M * M];
        for(int i = 0; i < M * M; ++i)
        {
           configuration[i] = 1;
        }
    }
    ~State()
    {
        free(configuration);
    }
    State(const State & state0)
    {
        M = state0.M;
        configuration = new int[M * M];
        for (int i = 0; i < M * M; ++i)
        {
            configuration[i] = state0.configuration[i];
        }
    }

    int get_spin(int i, int j)//获取spin的值
    {
        if(i >= M || j >= M || i < 0 || j < 0)
        {
            cout << "Warning! Index out of range!"<<endl;
            cout<<"("<<i<<","<<j<<")"<<endl;
            return 0;
        }
        return configuration[M * i + j];
    }

    double get_Hamilton()
    {
        double H = 0;
        for(int i = 0; i < M; ++i)
        {
            for(int j = 0; j < M; ++j)
            {
                H -= h * get_spin(i, j);
                if(i < M - 1)
                {
                    H -= J * (get_spin(i, j) * get_spin(i + 1, j));
                }
                if(j < M - 1)
                {
                    H -= J * (get_spin(i, j) * get_spin(i, j + 1));
                }
            }
        }
        return H;
    }

    double get_magnetization()
    {
        double mag = 0;
        for(int i = 0; i < M * M; i++)
            mag += configuration[i];
        return mag;
    }

    double get_correlation(int r)
    {
    	double gamma = 0.0;
    	if (r < 0 || r > M)
    		cout<<"Warning, r out of range."<<endl;

    	int rand_i = rand() % (M - r);
    	int rand_j = rand() % M;

    	gamma = get_spin(rand_i, rand_j) * get_spin(rand_i + r, rand_j);

    	return gamma;
    }


    void transition()
    {
        //Gibbs Sampling;
        //Proposal:
        int k = rand() % (M * M);//Gibbs采样中变化的分量；
        int ik = k / M;
        int jk = k % M;//ik和jk是该点对应的坐标；
        //cout<<"("<<ik<<","<<jk<<")"<<endl;

        //计算Hamildon量的变化；
        double H_old = 0;
        H_old -= h * get_spin(ik, jk);
        if(ik > 0)
            H_old -= J * (get_spin(ik, jk) * get_spin(ik - 1, jk));
        if(jk > 0)
            H_old -= J * (get_spin(ik, jk) * get_spin(ik, jk - 1));
        if(ik < M - 1)
            H_old -= J * (get_spin(ik, jk) * get_spin(ik + 1, jk));
        if(jk < M - 1)
            H_old -= J * (get_spin(ik, jk) * get_spin(ik, jk + 1));
        double H_new = -H_old;
        double dH = H_new - H_old;

        //Decision:
        if(dH < 0)
            configuration[k] = -configuration[k];
        else
        {
            double r = (rand() + 0.0) / RAND_MAX;
            if(r < exp(-beta * dH))
                configuration[k] = -configuration[k];
        }
//        double r = (rand() + 0.0) / RAND_MAX;
//        if(r < 1/(1 + exp(beta * dH)))
//            configuration[k] = -configuration[k];
//        return;
    }


};
int main() {
    srand(time(NULL));
    int M = 50;
    double EH = 0;
    double EH_square = 0;
    double Var_H = 0;
    double sum_H = 0;
    double sum_H_square = 0;
    double sum_mag = 0;
    double E_mag = 0;
    double dH = 0;

    double * Gamma = new double[M];
    double * EGamma = new double[M];

    int i = 0;
    double T = 0.1;

    ofstream fout("output.txt");

    while (T < 5)
    {
        State state = State(M, 1 / T, 1, 0);
        for (i = 0; i < 300000; ++i) {
            dH = state.get_Hamilton();
            sum_H += dH;
            sum_H_square += (dH * dH);

            sum_mag += state.get_magnetization();

            for(int r = 0; r < M; ++r)
            {
            	Gamma[r] += state.get_correlation(r);
            }

            state.transition();
        }

        EH = sum_H / (i + 1);
        EH_square = sum_H_square / (i + 1);
        Var_H = EH_square - EH * EH;

        E_mag = sum_mag / (i + 1);

        for (int r = 0; r < M; ++r)
        {
        	EGamma[r] = Gamma[r] / (i + 1);
        }
        cout <<"Temperature:"<<T<<endl;
        fout << T << " ";
        cout <<"Energy:"<< EH / (M * M) << endl;
        fout << EH / (M * M) << " ";
        cout << "Specific Heat:"<< Var_H / (T * T) / (M * M) << endl;
        fout << Var_H / (T * T) / (M * M) << " ";
        cout << "Magnetization:"<< E_mag / (M * M) << endl;
        fout << E_mag / (M * M) << " ";
        for (int r = 0; r < M; ++r)
        {
        	cout << "Gamma("<<r<<"):" << EGamma[r] <<endl;
        	fout << EGamma[r] << " ";
        }
        fout<<endl;

        sum_H = 0;
        sum_H_square = 0;
        sum_mag = 0;
        for (int r = 0; r < M; ++r)
        {
        	Gamma[r] = 0;
        	EGamma[r] = 0;
        }
        if (T >= 2 && T <= 2.495)
            T += 0.01;
        else
            T += 0.1;
    }

    free(Gamma);
    free(EGamma);
    fout.close();
    return 0;
}