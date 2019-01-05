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

    void print_configuration(ofstream & fout)
    {
        for (int i = 0; i < M; ++i)
        {
            for (int j = 0; j < M; ++j)
            {
                fout << get_spin(i, j)<<" ";
            }
            fout<<endl;
        }
    }


};
int main() {
    srand(time(NULL));
    int M = 50;
    double sum_mag = 0;
    double E_mag = 0;

    int i = 0;
    double h = 0;
    double Ts[6] = {0.1, 1, 2, 2.3, 3, 4};
    ofstream fout1("output1.txt");
    ofstream fout2("output2.txt");
    for (int j = 0; j < 6; ++j)
    {
        while (h < 3)
        {
            State state = State(M, 1 / Ts[j], 1, h);
            for (i = 0; i < 300000; ++i) {

                sum_mag += state.get_magnetization();

                state.transition();
            }

            E_mag = sum_mag / (i + 1);

            cout <<"External Magnetization:"<<h<<endl;
            fout1 << h << " ";
            cout << "Magnetization:"<< E_mag / (M * M) << endl;
            fout1 << E_mag / (M * M) << " ";

            sum_mag = 0;

            if(h <= 0.05)
                state.print_configuration(fout2);
            h += 0.1;
        }
        h = 0;
        fout1<<endl;
    }
    fout1.close();
    fout2.close();
    return 0;
}