#include <iostream>
#include <cmath>
#include <Eigen/Core>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

//
// 内点法による二次計画問題(Quadratic Programming: QP)の最適化
// QP:
// min 1/2 x^T Q x + c^T x
// s.t. A x = 0
//        x >= 0
//
void solveQPbyIPM(const MatrixXd& Q, const MatrixXd& A, const VectorXd& b, const VectorXd& c, VectorXd& x, VectorXd& nu) {

    int n = x.size(); // 変数の数
    int m = b.size(); // 制約の数
    double mu = 1; // ログバリアの強さを制御するパラメータ
    double beta = 0.1; // muの減衰率

    int maxOuterIter = 10; // muを減らす外側ループの上限
    int maxInnerIter = 100; // ニュートン法の繰り返し上限
    double tol = 1e-6; // 停止基準

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

            double f = (x.dot(Q*x))/2 + c.dot(x);
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
            // xの逆数を要素に持つベクトルd
            VectorXd d(x.size());
            for(int i=0; i<d.size(); i++)
                d(i) = 1 / x(i);
            // ラグランジュ関数のx微分
            Fy1 = 1/2 * (Q+Q.transpose())*x + c - mu*d - A.transpose()*nu;
            // 等式制約
            Fy2 = b - A*x;

            VectorXd Fy(Fy1.size()+Fy2.size());
            Fy << Fy1, Fy2; // 二つのベクトルを縦に結合

            // -- Fyのノルムがtol以下ならばinner loopを終了するコードを作成 --
            if(Fy.norm()<tol) { // 条件部を作成
                break;
            }



            // -- 1/x_i^2を対角に持つ行列を定義するコードを作成 --
            // n x n行列
            MatrixXd XXi = MatrixXd::Zero(n,n);
            for(int i=0; i<x.size(); i++)
                XXi(i, i) = 1 / (x(i)*x(i));

            // -- ヤコビ行列を計算するコードを作成 --
            // (n+m) x (n+m)行列になる
            // 右下はmxmの零行列
            MatrixXd J = MatrixXd::Zero(n+m, n+m);
            J << 1/2 * (Q+Q.transpose()) + mu*XXi,    -A.transpose(),
                    -A,                               MatrixXd::Zero(m,m);

            // -- ニュートン法の更新方向を計算するコードを作成 --
            VectorXd DeltaY(J.rows()); // xと\nuの更新方向（Jの行と同じ長さ）
            // ヤコビ行列の逆行列を用いて得られる
            DeltaY = -J.inverse()*Fy;

            // -- ステップ幅調整（fraction to the boundary ruleを参照）のコードを作成 --
            double tau = 0.995;
            double alpha=1;
            // i=0からnまで繰り返し
            for(int i=0; i<x.size(); i++) {
                double newAlpha = - (x(i)/DeltaY(i)) * tau;
                // newAlphaが条件を満たすときalphaを更新
                if(newAlpha>0 && alpha>newAlpha)
                    alpha = newAlpha;
            }

            // -- xと\nuをステップ幅\alphaで更新するコードを作成 --
            // DeltaYの上部がx,下部がnuなのでhead()とtail()でアクセス
            x  = x  + alpha*DeltaY.head(x.size());
            nu = nu + alpha*DeltaY.tail(nu.size());

        } // end of inner iteration

        // 変化量が小さくなったら終了
        if((x - xOld).squaredNorm() < tol) {
            break;
        }

        // 内側ループ（ニュートン法）が収束した段階での解
        cout << "Inner iteration finished" << endl; // rm
        cout << "x = " << x.transpose() << endl;

        // バリアパラメータの更新
        mu = beta*mu;

    } // end of outer iteration

}

int main() {

    int n = 6; // 変数の数 x: n次元ベクトル
    int m = 4; // 等式制約の数 A: (n, m) 次元, b: m次元ベクトル

    // 係数を設定する行列とベクトル
    MatrixXd Q(n,n);
    VectorXd c(n);
    MatrixXd A(m,n);
    VectorXd b(m);

    Q <<
            2, 1, 0, 0, 0, 0,
            1, 2, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0;
    A <<
            1,-1, 1, 0, 0, 0,
            1, 3, 0, 1, 0, 0,
            3, 1, 0, 0, 1, 0,
            -3, 1, 0, 0, 0, 1;
    b << 3, 5, 8, 10;
    c << -6,-6, 0, 0, 0, 0;

    VectorXd x(n);
    VectorXd nu(m);
     x << 1, 1, 1, 1, 1, 1; // 等式制約と非負制約を満たす適当な初期値を設定
    nu << 0, 0, 0, 0; // 双対変数の初期値

    solveQPbyIPM(Q, A, b, c, x, nu);

    cout << "==================================================" << endl;
    cout << "Optimal solution:" << endl;
    cout << " x = " << x.transpose() << endl;

    return 0;
}
