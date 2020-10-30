# Internal Point Method

内点法(Internal Point Method)で最適化問題を解いてみた．

### 線形計画問題

`asg01/main.cpp`では以下の標準形の線形計画問題を内点法によって解く．

![
\begin{align*}
\min_{ x} \  &  c^\top x\\
 \text{subject to} \ & A x = b \\
 & x \geq 0
\end{align*}](https://render.githubusercontent.com/render/math?math=%5Clarge+%5Cdisplaystyle+%0A%5Cbegin%7Balign%2A%7D%0A%5Cmin_%7B+x%7D+%5C++%26++c%5E%5Ctop+x%5C%5C%0A+%5Ctext%7Bsubject+to%7D+%5C+%26+A+x+%3D+b+%5C%5C%0A+%26+x+%5Cgeq+0%0A%5Cend%7Balign%2A%7D)

例えば，以下の問題

![
\begin{align*}
\min_{x_1, x_2} \ & -5 x_1 - 4 x_2 \\
 \text{subject to} \ 
 & 5 x_1 + 2 x_2 \leq 30 \\
 &   x_1 + 2 x_2 \leq 14 \\
 & 5 x_1 - 4 x_2 \leq 15 \\
 & 5 x_1 - 2 x_2 \leq 20 \\
 &   x_1 \geq 0, \ x_2 \geq 0
\end{align*}](https://render.githubusercontent.com/render/math?math=%5Clarge+%5Cdisplaystyle+%0A%5Cbegin%7Balign%2A%7D%0A%5Cmin_%7Bx_1%2C+x_2%7D+%5C+%26+-5+x_1+-+4+x_2+%5C%5C%0A+%5Ctext%7Bsubject+to%7D+%5C+%0A+%26+5+x_1+%2B+2+x_2+%5Cleq+30+%5C%5C%0A+%26+++x_1+%2B+2+x_2+%5Cleq+14+%5C%5C%0A+%26+5+x_1+-+4+x_2+%5Cleq+15+%5C%5C%0A+%26+5+x_1+-+2+x_2+%5Cleq+20+%5C%5C%0A+%26+++x_1+%5Cgeq+0%2C+%5C+x_2+%5Cgeq+0%0A%5Cend%7Balign%2A%7D)

は，以下のように変換することでこのプログラムで解くことができる．

![
\begin{align*}
\min_{x_1, x_2} \ & -5 x_1 - 4 x_2 \\
 \text{subject to} \ 
 & 5 x_1 + 2 x_2 + x_3 = 30 \\
 &   x_1 + 2 x_2 + x_4 = 14 \\
 & 5 x_1 - 4 x_2 + x_5 = 15 \\
 & 5 x_1 - 2 x_2 + x_6 = 20 \\
 &   x_1, x_2, x_3, x_4, x_5, x_6 \geq 0
\end{align*}](https://render.githubusercontent.com/render/math?math=%5Clarge+%5Cdisplaystyle+%0A%5Cbegin%7Balign%2A%7D%0A%5Cmin_%7Bx_1%2C+x_2%7D+%5C+%26+-5+x_1+-+4+x_2+%5C%5C%0A+%5Ctext%7Bsubject+to%7D+%5C+%0A+%26+5+x_1+%2B+2+x_2+%2B+x_3+%3D+30+%5C%5C%0A+%26+++x_1+%2B+2+x_2+%2B+x_4+%3D+14+%5C%5C%0A+%26+5+x_1+-+4+x_2+%2B+x_5+%3D+15+%5C%5C%0A+%26+5+x_1+-+2+x_2+%2B+x_6+%3D+20+%5C%5C%0A+%26+++x_1%2C+x_2%2C+x_3%2C+x_4%2C+x_5%2C+x_6+%5Cgeq+0%0A%5Cend%7Balign%2A%7D)

以下は実行例．

```
==================================================
mu = 1
--------------------------------------------------
Innter iteration: 0
f(x) = -9
f(x) - mu \sum log x = -9
A x - b = -22 -10 -13 -16
x = 1 1 1 1 1 1
nu = 0 0 0 0
--------------------------------------------------
Innter iteration: 1
f(x) = -35.2761
f(x) - mu \sum log x = -30.4239
A x - b = 3.55271e-15           0           0           0
x =  4.75118  2.88004 0.484023  3.48874  2.76427  2.00419
nu =   -1.51598    1.48874   0.764274 0.00419068
--------------------------------------------------
(略)
==================================================
Optimal solution:
 x =           4           5 1.33334e-05 7.99996e-06          15          10
```

`mu`は小さいほど精度が良くなる．

### 2次計画問題

`asg02/main.cpp`では以下の標準形の2次計画問題を内点法によって解く．

![
\begin{align*}
 \min_{x} \ & \frac{1}{2} x^\top Q x + c^\top x \\
  \text{subject to} \ & A x = b \\
  & x \geq 0
\end{align*}](https://render.githubusercontent.com/render/math?math=%5Clarge+%5Cdisplaystyle+%0A%5Cbegin%7Balign%2A%7D%0A+%5Cmin_%7Bx%7D+%5C+%26+%5Cfrac%7B1%7D%7B2%7D+x%5E%5Ctop+Q+x+%2B+c%5E%5Ctop+x+%5C%5C%0A++%5Ctext%7Bsubject+to%7D+%5C+%26+A+x+%3D+b+%5C%5C%0A++%26+x+%5Cgeq+0%0A%5Cend%7Balign%2A%7D)

例えば，以下の問題

![
\begin{align*}
 \begin{split}  
\min_{x_1, x_2} \ &  \frac{1}{2} 
\begin{bmatrix}
  x_1 & x_2
 \end{bmatrix}
 \begin{bmatrix}
  2 & 1 \\
  1 & 2
 \end{bmatrix}
 \begin{bmatrix}
  x_1 \\
  x_2
 \end{bmatrix}+ 
 \begin{bmatrix}
  -6 & -6
 \end{bmatrix}
 \begin{bmatrix}
  x_1 \\
  x_2
 \end{bmatrix}
\\
 \text{subject to} \ 
  & x_1 - x_2 \leq 3 \\
  & x_1 + 3 x_2 \leq 5 \\
  & 3 x_1 + x_2 \leq 8 \\
  & -3 x_1 + x_2 \leq 10 
 \end{split}
\end{align*}](https://render.githubusercontent.com/render/math?math=%5Clarge+%5Cdisplaystyle+%0A%5Cbegin%7Balign%2A%7D%0A+%5Cbegin%7Bsplit%7D++%0A%5Cmin_%7Bx_1%2C+x_2%7D+%5C+%26++%5Cfrac%7B1%7D%7B2%7D+%0A%5Cbegin%7Bbmatrix%7D%0A++x_1+%26+x_2%0A+%5Cend%7Bbmatrix%7D%0A+%5Cbegin%7Bbmatrix%7D%0A++2+%26+1+%5C%5C%0A++1+%26+2%0A+%5Cend%7Bbmatrix%7D%0A+%5Cbegin%7Bbmatrix%7D%0A++x_1+%5C%5C%0A++x_2%0A+%5Cend%7Bbmatrix%7D%2B+%0A+%5Cbegin%7Bbmatrix%7D%0A++-6+%26+-6%0A+%5Cend%7Bbmatrix%7D%0A+%5Cbegin%7Bbmatrix%7D%0A++x_1+%5C%5C%0A++x_2%0A+%5Cend%7Bbmatrix%7D%0A%5C%5C%0A+%5Ctext%7Bsubject+to%7D+%5C+%0A++%26+x_1+-+x_2+%5Cleq+3+%5C%5C%0A++%26+x_1+%2B+3+x_2+%5Cleq+5+%5C%5C%0A++%26+3+x_1+%2B+x_2+%5Cleq+8+%5C%5C%0A++%26+-3+x_1+%2B+x_2+%5Cleq+10+%0A+%5Cend%7Bsplit%7D%0A%5Cend%7Balign%2A%7D)

は，以下のように変換することでこのプログラムで解くことができる．

![\begin{align*}
\begin{split}  
 \min_{x_1, x_2} \ & \frac{1}{2} 
\begin{bmatrix}
  x_1 & x_2
 \end{bmatrix}
 \begin{bmatrix}
  2 & 1 \\
  1 & 2
 \end{bmatrix}
 \begin{bmatrix}
  x_1 \\
  x_2
 \end{bmatrix} +
\begin{bmatrix}
  -6 & -6
 \end{bmatrix}
 \begin{bmatrix}
  x_1 \\
  x_2
 \end{bmatrix}
 \\
 \text{subject to} \ 
  & x_1 - x_2 + x_3 = 3 \\
  & x_1 + 3 x_2 + x_4 = 5 \\
  & 3 x_1 + x_2 + x_5 = 8 \\
  & -3 x_1 + x_2 + x_6 = 10 \\
  & x_3, x_4, x_5, x_6 \geq 0
 \end{split}
\end{align*}](https://render.githubusercontent.com/render/math?math=%5Clarge+%5Cdisplaystyle+%5Cbegin%7Balign%2A%7D%0A%5Cbegin%7Bsplit%7D++%0A+%5Cmin_%7Bx_1%2C+x_2%7D+%5C+%26+%5Cfrac%7B1%7D%7B2%7D+%0A%5Cbegin%7Bbmatrix%7D%0A++x_1+%26+x_2%0A+%5Cend%7Bbmatrix%7D%0A+%5Cbegin%7Bbmatrix%7D%0A++2+%26+1+%5C%5C%0A++1+%26+2%0A+%5Cend%7Bbmatrix%7D%0A+%5Cbegin%7Bbmatrix%7D%0A++x_1+%5C%5C%0A++x_2%0A+%5Cend%7Bbmatrix%7D+%2B%0A%5Cbegin%7Bbmatrix%7D%0A++-6+%26+-6%0A+%5Cend%7Bbmatrix%7D%0A+%5Cbegin%7Bbmatrix%7D%0A++x_1+%5C%5C%0A++x_2%0A+%5Cend%7Bbmatrix%7D%0A+%5C%5C%0A+%5Ctext%7Bsubject+to%7D+%5C+%0A++%26+x_1+-+x_2+%2B+x_3+%3D+3+%5C%5C%0A++%26+x_1+%2B+3+x_2+%2B+x_4+%3D+5+%5C%5C%0A++%26+3+x_1+%2B+x_2+%2B+x_5+%3D+8+%5C%5C%0A++%26+-3+x_1+%2B+x_2+%2B+x_6+%3D+10+%5C%5C%0A++%26+x_3%2C+x_4%2C+x_5%2C+x_6+%5Cgeq+0%0A+%5Cend%7Bsplit%7D%0A%5Cend%7Balign%2A%7D)

以下は実行例．

```
==================================================
mu = 1
--------------------------------------------------
Innter iteration: 0
f(x) = -9
f(x) - mu \sum log x = -9
A x - b =  -2   0  -3 -11
x = 1 1 1 1 1 1
nu = 0 0 0 0
--------------------------------------------------
Innter iteration: 1
f(x) = -9.21562
f(x) - mu \sum log x = -11.5215
A x - b =  -1.3275        0 -1.99125 -7.30125
x = 0.68625 1.43625  2.4225   0.005 2.51375 3.32125
nu =  1.08625 -1.33125   1.1775    1.985
--------------------------------------------------
(略)
==================================================
Optimal solution:
 x =     2.37498    0.874983         1.5 6.66637e-05 6.66677e-05       16.25
```

上と同様に，`mu`が小さくなるほど精度が良くなる．