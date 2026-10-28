#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <vector>
using namespace std;


const int N = 8297;
const double M = 1.0 / N;
double G[N + 1][N + 1], X[2][N + 1], p = 0.85;
int C[N + 1], R[N + 1], Cnt;
bool Ch = 0;


double Auth[N + 1], Hub[N + 1];
bool Er;

// 8297
// 103689

void CalG()
{
    for (int i = 1; i <= N; i++)
        for (int j = 1; j <= N; j++)
            G[i][j] = p * G[i][j] + (1 - p) * M;
}

void Mul()
{
    for (int i = 1; i <= N; i++)
    {
        double Result = 0.0;
        int Val = Cnt & 1;
        for (int j = 1; j <= N; j++) Result += G[i][j] * X[Val][j];
        X[1 - Val][i] = Result;
    }
}

double Error(double* A, double* B)
{
    double Result = 0.0;
    for (int i = 1; i <= N; i++) Result += (A[i] - B[i]) * (A[i] - B[i]);
    return sqrt(Result);
}

void Normalize()
{
    double Sum1 = 0.0;
    double Sum2 = 0.0;
    for (int i = 1; i <= N; i++)
    {
        Sum1 += Auth[i] * Auth[i];
        Sum2 += Hub[i] * Hub[i];
    }

    double Norm1 = sqrt(Sum1);
    double Norm2 = sqrt(Sum2);
    for (int i = 1; i <= N; i++)
    {
        if (Norm1 > 0) Auth[i] /= Norm1;
        if (Norm2 > 0) Hub[i] /= Norm2;
    }
}

void Export()
{
    ofstream Out("edges.txt");
    for (int i = 1; i <= N; i++)
        for (int j = 1; j <= N; j++)
            if (C[j] > 0 && G[i][j] > (1.0 - p) * M + 1e-15) Out << j << ' ' << i << '\n';
    Out.close();
}

double Angl(double* A, double* B)
{
    double Dot = 0.0, Sum1 = 0.0, Sum2 = 0.0;
    for (int i = 1; i <= N; i++)
    {
        Dot += A[i] * B[i];
        Sum1 += A[i] * A[i];
        Sum2 += B[i] * B[i];
    }

    Sum1 = sqrt(Sum1);
    Sum2 = sqrt(Sum2);

    double Cs = Dot / (Sum1 * Sum2);
    if (Cs > 1) Cs = 1;
    if (Cs < -1) Cs = -1;

    double rads = acos(Cs);

    return rads * (180.0 / M_PI);
}

int main()
{
    //----------------------------------------------------- INPUT
    ifstream Fl("wiki-Vote.txt");
    string Line;
    while (getline(Fl, Line))
    {
        stringstream Str(Line);
        int u, v;
        Str >> u >> v;
        G[v][u] += 1.0;
        C[u] += 1;
        R[v] += 1;
    }

    cout << "Memory: " << N * N * 8.0 / (1024 * 1024) << " MG" << endl;
    cout << "Sparsity: " << 1.0 - 103689.0 / (N * N) << endl;

    bool Sym = 1;
    for (int i = 1; i <= N; i++)
        for (int j = 1; j <= N; j++)
            if (G[i][j] != G[j][i]) Sym = 0;
    
    if (!Sym) cout << "Not Symmetrical." << endl;
    else cout << "Symmetrical." << endl; 

    cout << endl;

    //----------------------------------------------------- Page Ranking

    for (int i = 1; i <= N; i++)
    {
        for (int j = 1; j <= N; j++)
        {
            if (C[j] == 0)
            {
                G[i][j] = M;
                continue;
            }

            G[i][j] /= C[j];
        }
    }
    
    CalG();

    for (int i = 1; i <= N; i++) X[0][i] = M;
    while (!Ch)
    {
        Mul();
        Ch = (Error(X[0], X[1]) < 1.0 / 1000000);
        Cnt++;
    }

    vector<pair<double , int> > Vec;
    for (int i = 1; i <= N; i++) Vec.push_back({X[Cnt & 1][i], i});
    sort(Vec.rbegin(), Vec.rend());
    cout << "Page Ranking: " << endl;
    for (int i = 0; i < 5; i++) cout << i + 1 << ": " << Vec[i].second << " (Score: " << Vec[i].first << ")" << endl;

    //----------------------------------------------------- HITS

    for (int i = 1; i <= N; i++)
    {
        Auth[i] = 1.0;
        Hub[i] = 1.0;
    }


    while (!Er)
    {
        double New_Auth[N + 1] = {0};
        double New_Hub[N + 1] = {0};
        double Old_Auth[N + 1] = {0};
        double Old_Hub[N + 1] = {0};

        for (int i = 1; i <= N; i++)
            for (int j = 1; j <= N; j++)
                if (C[j] != 0 && G[i][j] > (1.0 - p) * M + 1e-15) New_Auth[i] += Hub[j];
        
        for (int i = 1; i <= N; i++)
            for (int j = 1; j <= N; j++)
                if (C[i] != 0 && G[j][i] > (1.0 - p) * M + 1e-15) New_Hub[i] += New_Auth[j];
        
        for (int i = 1; i <= N; i++)
        {
            Old_Auth[i] = Auth[i];
            Old_Hub[i] = Hub[i];
            Auth[i] = New_Auth[i];
            Hub[i] = New_Hub[i];
        }

        Normalize();

        Er = (max(Error(Auth, Old_Auth), Error(Hub, Old_Hub)) < 1.0 / 10000000);
    }

    vector<pair<double , int> > AuL, HuL;
    for (int i = 1; i <= N; i++)
    {
        AuL.push_back({Auth[i], i});
        HuL.push_back({Hub[i], i});
    }

    sort(AuL.rbegin(), AuL.rend());
    sort(HuL.rbegin(), HuL.rend());

    cout << endl;

    cout << "Authority Ranking: " << endl;
    for (int i = 0; i < 5; i++)
        cout << i + 1 << ": " << AuL[i].second << " (Score: " << AuL[i].first << ") " << endl;

    cout << endl;

    cout << "Hub Ranking: " << endl;
    for (int i = 0; i < 5; i++)
        cout << i + 1 << ": " << HuL[i].second << " (Score: " << HuL[i].first << ") " << endl;
    
    //----------------------------------------------------- Comparing

    vector<pair<int, int> > Deg;
    for (int i = 1; i <= N; i++)
        Deg.push_back({R[i], i});
    sort(Deg.rbegin(), Deg.rend());

    cout << endl;

    cout << "Degree Ranking: " << endl;
    for (int i = 0; i < 5; i++)
        cout << i + 1 << ": " << Deg[i].second << " (Score: " << Deg[i].first << ") " << endl;

    Export();

    system("python Help.py");

    ifstream In("svd_results.txt");
    vector<pair<double, int> > AnsAu, AnsHu;
    double Atu[N + 1], Hb[N + 1];
    string Ln;
    getline(In, Ln);
    stringstream S1(Ln);
    for (int i = 1; i <= N; i++)
    {
        S1 >> Atu[i];
        AnsAu.push_back({abs(Atu[i]), i});
    }

    getline(In, Ln);
    stringstream S2(Ln);
    for (int i = 1; i <= N; i++)
    {
        S2 >> Hb[i];
        AnsHu.push_back({abs(Hb[i]), i});
    }

    getline(In, Ln);
    stringstream S3(Ln);
    double Value;
    S3 >> Value;

    sort(AnsAu.rbegin(), AnsAu.rend());

    sort(AnsHu.rbegin(), AnsHu.rend());

    cout << endl;

    cout << "Greatest singular value of Adj Matrix: " << Value << endl;
    cout << "Greatest eigenvalue of Adj Matrix: " << Value * Value << endl;

    cout << endl;

    cout << "Real Authority TOP Five: " << endl;
    for (int i = 0; i < 5; i++)
        cout << i + 1 << ": " << AnsAu[i].second << " (Score: " << AnsAu[i].first << ") " << endl;

    cout << endl;

    cout << "Real Hubs TOP Five: " << endl;
    for (int i = 0; i < 5; i++)
        cout << i + 1 << ": " << AnsHu[i].second << " (Score: " << AnsHu[i].first << ") " << endl;

    cout << endl;

    cout << "Angel between U1 and our Authority Ranking: " << Angl(Atu, Auth) << endl;
    cout << "Angel between V1 and our Hub Ranking: " << Angl(Hb, Hub) << endl;

    cout << endl;
    return 0;
}