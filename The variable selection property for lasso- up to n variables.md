# The variable selection property for lasso-> up to n variables

For the lasso form penalty problem 
$$
\hat { \beta } \in \underset { \beta \in \mathbb { R } ^ { p } } { \operatorname { argmin } } \frac { 1 } { 2 } \| y - X \beta \| _ { 2 } ^ { 2 } + \lambda \| \beta \| _ { 1 }
$$
the elastic net paper(2004) said that the lasso can only select up to n parameters. Here are two way to explain this.

1. The optimization way. Quote from [stackexchange][https://stats.stackexchange.com/questions/38299/if-p-n-the-lasso-selects-at-most-n-variables]

   The lasso penalty is a optimization problem which fits the KKT condition, so we have the KKT conditon as
   $$
   \left. \begin{array} { c } { X ^ { T } ( y - X \hat { \beta } ) = \lambda \gamma } \\ { \gamma _ { i } \in \left\{ \begin{array} { l l } { \left\{ \operatorname { sign } \left( \hat { \beta } _ { i } \right) \right\} } & { \text { if } \hat { \beta } _ { i } \neq 0 } \\ { [ - 1,1 ] } & { \text { if } \hat { \beta } _ { i } = 0 } \end{array} \right. , \text { for } i = 1 , \ldots p } \end{array} \right.
   $$
   We use the main of the KKT conditon that $|X^T(y-X\hat\beta)|=\lambda$, in the abosolute value case we can omit the sub-gradient.

   This formula is a linear system with $\hat\beta$. Since the rank of the matrix $X$ is up to n, so this linear system with unique result only with n line of $\beta$ . Otherwise it will have infinity solutions. ​	

   In another words, this type of optimization problem only solvable when active set $\Epsilon$ is not larger than n predictors.	

   So in such point of view, the key thing is the rank about the hat matrix. In OLS is $X^TX$. 													

2. The approximation to OLS view of the lasso problem.(by chunyu)

   When lambda is small and close to 0, the lasso estimator is come close to OLS which is solvable only up to n estimators.

   The problem is, ridge regression is also close to OLS when $\lambda$ comes to 0. Why this unsolvable situation do not happen to ridgeregression ?

   The answer is, ridge regression always keep all the predictor in the model, however, for lasso, the number of predictors be selected is monotonic decreaasing via $\lambda$. The number of predictor is maximized when $\lambda$ close to 0, that is, the OLS case, number of predictors equals to n.

3. The LARS point of view about this

   Tbc..

In the mean time, the ridge-regression form is as following
$$
\hat\beta=(X^TX+\lambda I)^{-1}X^Ty
$$
The rank of hat matrix $X^TX+\lambda I$ can equal to n, so it is solvable.

And for the KKT condition case, maybe because it is strictly convex, so the form of $\gamma$ and the equation we use could change to
$$
\begin{align}
X^T(y-X\hat\beta)=2\lambda\hat\beta\\
(X^TX+\lambda I)\hat\beta=X^Ty
\end{align}
$$


Then the rank of the linear equations is not up to n but up to p. So for elastic net is same, but combine these two matrix together.









