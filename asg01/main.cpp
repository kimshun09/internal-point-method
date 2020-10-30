#include <iostream>
#include <cmath>
#include <Eigen/Core>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

//
// 内点法(Interior Point Method: IPM)による線形計画問題(Linear Programming: LP)の最適化
// LP (標準形):
// min c^T x
// s.t. A x = b
//        x >= 0
//

void solveLPbyIPM(const MatrixXd& A, const VectorXd& b, const VectorXd& c, VectorXd& x, VectorXd& nu) {

    int n = x.size(); // 変数の数
    int m = b.size(); // 制約の数
    double mu = 1; // ログバリアの強さを制御するパラメータ
    double beta = 0.1; // muの減衰率

    int maxOuterIter = 10; // muを減らす外側ループの上限
    int maxInnerIter = 100; // ニュートン法の繰り返し上限
    double tol = 1e-6; // 停止基準 (tolerance)

    //
    // Outer loop
    //
    for(int outerIter = 0; outerIter < maxOuterIter; ++outerIter) {
        cout << "==================================================" << endl;
        cout << "mu = " << mu << endl;
        VectorXd xOld = x;

        //
        // Innter loop:
        // ある\muに関するログバリア最適化問題の
        // KKT条件の方程式をNewton法で解く
        //
        for(int innerIter = 0; innerIter < maxInnerIter; ++innerIter) {

            double f = c.dot(x);
            double logX = x.array().log().sum();

            cout << "--------------------------------------------------" << endl;
            cout << "Innter iteration: " << innerIter << endl;
            cout << "f(x) = " << f << endl; // 目的関数値
            cout << "f(x) - mu \\sum log x = " << f + mu * logX << endl; // 目的関数値 + logペナルティ
            cout << "A x - b = " << (A * x - b).transpose() << endl; // 制約条件
            cout << "x = " << x.transpose() << endl; // 主変数
            cout << "nu = " << nu.transpose() << endl; // 双対変数

            //
            // KKT条件の方程式を計算
            //
            // F(y) = [ Fy1 ]
            //        [ Fy2 ]
            //
            VectorXd Fy1(n); // ラグランジュ関数のxでの微分: d L(x,\nu) / d x
            VectorXd Fy2(m); // 等式制約

            // -- Fy1とFy2を計算するコードを作成 --
            // 要素がxの逆数になっているベクトルd
            VectorXd d(x.size());
            for(int i=0; i<d.size(); i++)
                d(i) = 1 / x(i);
            // ラグランジュ関数のx微分
            Fy1 = c - mu*d - A.transpose()*nu;
            // 等式制約
            Fy2 = b - A*x;


            VectorXd Fy(Fy1.size()+Fy2.size());
            Fy << Fy1, Fy2; // 二つのベクトルを縦に結合

            // -- Fyのノルムが変数tol以下ならばinner loopを終了するコードを作成 --
            if(Fy.norm()<tol) { // 条件部を作成
                break;
            }

            // -- 1/x_i^2を対角に持つ行列を定義するコードを作成 --
            // n x n行列になる
            MatrixXd XXi = MatrixXd::Zero(n,n);
            for(int i=0; i<x.size(); i++)
                XXi(i, i) = 1 / (x(i)*x(i));

            // -- ヤコビ行列を計算するコードを作成 --
            // (n+m) x (n+m)行列になる
            // 右下は0ではなく零行列
            MatrixXd J = MatrixXd::Zero(n+m, n+m);
            J << mu*XXi, -A.transpose(),
                 -A,     MatrixXd::Zero(m, m);

            // -- ニュートン法の更新方向を計算するコードを作成 --
            VectorXd DeltaY(J.rows()); // xと\nuの更新方向（Jの行と同じ長さ）
            // ヤコビ行列の逆行列をもちいる
            DeltaY = -J.inverse() * Fy;

            // -- ステップ幅調整（fraction to the boundary ruleを参照）のコードを作成 --
            double tau = 0.995;
            double alpha=1;
            // i=0からnまで繰り返し
            for(int i=0; i<x.size(); i++) {
                double newAlpha = - (x(i)/DeltaY(i)) * tau;
                // newAlphaが条件を満たせばalphaを更新
                if(newAlpha>0 && alpha>newAlpha)
                    alpha = newAlpha;
            }

            // -- xと\nuをステップ幅\alphaで更新するコードを作成 --
            // DeltaYは上部にx,下部にyが格納されているので
            // head()またはtail()で抜き出せる
            x  = x  + alpha*DeltaY.head(x.size());
            nu = nu + alpha*DeltaY.tail(nu.size());

        } // end of inner iteration

        // 変化量が小さくなったら終了
        if((x - xOld).squaredNorm() < tol) {
            break;
        }

        // 内側ループ（ニュートン法）が終了した段階での解
        cout << "Inner iteration finished" << endl;
        cout << "x = " << x.transpose() << endl;

        // バリアパラメータの更新
        mu = beta*mu;

    } // end of outer iteration

}


int main() {

    int n = 6; // 変数の数 x_1, ..., x_n
    int m = 4; // 等式制約の数 h_1, ..., h_m

    //
    // 線形計画問題
    //
    // min_x c^T x
    // s.t.  A x = b
    //         x >= 0
    //
    VectorXd c(n);
    MatrixXd A(m,n);
    VectorXd b(m);

    A << // 例題の係数
            5, 2, 1, 0, 0, 0,
            1, 2, 0, 1, 0, 0,
            5,-4, 0, 0, 1, 0,
            5,-2, 0, 0, 0, 1;
    b << 30, 14, 15, 20;
    c << -5, -4, 0, 0, 0, 0;

    cout << "LP parameters:";
    cout << A << endl;
    cout << b << endl;
    cout << c << endl;

    VectorXd x(n);
    VectorXd nu(m);
    x << 1, 1, 1, 1, 1, 1; // 初期値
    nu << 0, 0, 0, 0;

    solveLPbyIPM(A, b, c, x, nu);

    cout << "==================================================" << endl;
    cout << "Optimal solution:" << endl;
    cout << " x = " << x.transpose() << endl;

    return 0;
}

