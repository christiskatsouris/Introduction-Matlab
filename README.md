# Introduction-Matlab
Teaching Page (Drafted July 2022).

TIP: Learn, Unlearn & Relearn.

- An interesting article from the World Economic Forum: [Why learning to unlearn prepares us for transformative change](https://www.weforum.org/agenda/2022/02/learning-to-unlearn-exponetially/). Similar principles apply when learning a new programming language or when you apply the "asymptotic theory laws" you learn in an Econometric Theory class in order to derive the limit theory of statistics and estimators. Revising what you know, and building on existing acquisition of knowledge is a good practise for programming as well. 

# 1. Basics of Matlab

To begin with, Matlab is an Object-Oriented Programming Language which has faster computational capabilities in comparison to other Statistical Software. The best way to master a programming language such as Matlab or C++ is to consider reviewing some key concepts from about 50 years ago. Two seminal references are: 

- Tenenbaum, A. M., & Augenstein, M. J. (1986). Data structures using Pascal. Prentice-Hall, Inc.
- Standish, T. A. (1980). Data structure techniques. Addison-Wesley Longman Publishing Co., Inc.

In summary the main idea behind this is to aim to produce clean coding procedures and reproducible Matlab Scripts. This can be especially helpful when reviewing your code as well as when developing MATLAB code to implement simuation studies (such as Monte Carlo simulation experiments, bootstrap resampling methods etc) or other statistical algorithms with high computational/execution time. Furthermore, this can be extremely helpful when colloborating on research projects with others.  

An Introduction to Matlab Programming can be found [here](https://uk.mathworks.com/academia/courseware/introduction-to-matlab.html).  


# 1.1. Loops and Logical Branching (Syntax)


### The for loop

A simple form of such a loop is 

```Matlab

for index = 1:n
  statements 
end

```

### The if-elseif-else statement

A simple form of the if statement is 

```Matlab

if (condition)
  statements 
end

```

A more general form is

```Matlab

if (condition1)
  statementsA
elseif (condition2)
  statementsB
elseif (condition3)
  ...
else 
 statementsE

```

### The while loop

```Matlab

while (condition)
  statements 
end

```

## Example 1

Suppose we have invested some money in  fund which pays 5% (compound) interest per year, and we would like to know how long it takes for the value of the investment to double. 

Remark: In practise this example demonstrates the main difference of the functionality of the 'for' loop versus the 'while' loop for looping purposes. In particular for the above calculation we need to obtain a statement of the account for each year until the balance is doubled. Therefore, we cannot use a 'for' loop in this case, because we do not know beforehand how long this will take, so we cannot assign a value for the number of iterations on entering the loop. Thus, we should use a 'while' loop. 

```Matlab

% Example: Financial Calculations 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%script tutorial.m

format bank 
invest = ('type initial investment:  ')

r    = 0.05;
bal  = invest;
year = 0;
disp('         Year        Balance')

while (bal < 2 * invest)

  bal   = bal + 2*ball;
  year  = year + 1;
  disp([year,bal])
  
end  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

```

# 1.2. Using Matrices and Vectors

In cases that the coding procedure has several replications, for instance in Monte Carlo simulation studies or when applying bootstrap resampling methods, the initial step is to initialize the vector and matrices to store the results of the simulation design under examination. 

```Matlab

H0 = [0, 0.1, 0.5;  % The null configurations
      0, 0.5, 2;
      0, 0.5, 1];

H1 = [a;b;c];    

DGP = [H0, H1];     % The full configurations
D = size(DGP,2);
a = DGP(1,:);
b = DGP(2,:);
c = DGP(3,:);

```

# 2. Numerical Analysis/Optimization Examples

### Cholesky Decomposition 

```Matlab

>> A = [2, -1, 0; -1,2,-1; 0, -1, 2 ]

A =

     2    -1     0
    -1     2    -1
     0    -1     2

>> chol(A)

ans =

    1.4142   -0.7071         0
         0    1.2247   -0.8165
         0         0    1.1547

```

### Gauss Approximation Method

```Matlab

function[x] = gaussel(A,b)

  N = max(size(A));
  
  for j = 2:N
    for i = j:N
        m = A(i,j-1) / A(j-1,j-1);
        A(i,:) = A(i,:) - A(j-1,:)*m;
        b(i)   = b(i) - m*b(j-1);
    end    
  end
  
  x    = zeros(N,1);
  x(N) = b(N) / A(N,N);
  for j = N-1:-1:1
     x(j) = ( b(j) - A(j,j+1:N)*x(j+1:N) ) / A(j,j);
  end
  
>> gassel(A,b)

% Reference: Γεωργίου and Ξενοφώντος (2007).

```


# 3. Programming in Matlab using functions 

Constructing a function (i.e., a small program) in Matlab allow us to organise our coding procedure more efficiently. Furthermore, 'calling' functions in the main workflow is a good programming practice and permits to break-down a procedure into various estimation steps. 

##  3.1. Syntax

For the Syntax see: https://uk.mathworks.com/help/matlab/ref/function.html 

```Matlab

function [y1,...,yN] = myfun(x1,...,xM)

% function [y1,...,yN] = myfun(x1,...,xM) declares a function named myfun 
% that  accepts inputs x1,...,xM and returns outputs y1,...,yN. 
% This declaration statement must be the first executable line of the function. 
% Valid function names begin with an alphabetic character and can contain letters, numbers, or underscores.

%%%%%%%%%%%%%%%%%%%%%%%%
%% Example 
%%%%%%%%%%%%%%%%%%%%%%%%

function ave = average(x)
    ave = sum(x(:))/numel(x); 
end

z = 1:99;
ave = average(z)

ave =
    50

```

## Example 2

A helpful programming practice to identify mistakes with either the coding procedure or with the use of parameters which are not allowed by definition when executing the program is to use output messages such as below. 

```Matlab

 if mod( i_sim,100 ) == 0
        disp([label ', pseudo, i_sim = ' num2str(i_sim)])
    end

```

## Example 3

Consider constructing a function that obtains the parameter estimates from a quantile regression model.    

```Matlab

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function beta = quantile_regression(y,x,k)

    [xn,xm] = size(x);

    x = [ones(xn,1) x];
    xm = xm + 1;
    xs = x;

    beta = ones(xm,1);

    diff = 1;
    iter = 0;

    while ((diff > 1e-6) && (iter < 1000))
        xst = xs';
        beta_0 = beta;

        beta = ((xst * x) \ xst) * y;

        rsd = y - (x * beta);
        rsd(abs(rsd) < 0.000001) = 0.000001;
        rsd(rsd < 0) = k * rsd(rsd < 0);
        rsd(rsd > 0) = (1 - k) * rsd(rsd > 0);
        rsd = abs(rsd);

        z = zeros(xn,xm);

        for i = 1:xm 
            z(:,i) = x(:,i) ./ rsd;
        end

        xs = z;
        beta_1 = beta;
        
        diff = max(abs(beta_1 - beta_0));
        iter = iter + 1;
    end
end % end-of-function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% References: Matlab Package 'SystemicRisk' and R package 'quantreg'. 

```

##  3.2. Simulation Experiments Examples

```Matlab

% Choose substreams for reproducibility
set(stream,'Substream',r);
RandStream.setGlobalStream(stream);

% Simulate a regressor Z such that 
Z = 1 + beta*normcdf( randn(N,1) );  

% Simulate an error sequence
U = randn(N,1); 

```

## Example 4

Suppose the above regressor and error term are the components of an econometric model of interest or in other words the data generating process (DGP). Now assume that we are interested to evaluate the empirical size of a test statistic Tn based on the underline DGP. In addition, we have 3 different experimental designs and in each case we shall compute and store the p-value of the test statistic that corresponds to the number of replications used in the simulation study. Following good programming and good parallelism practices in Matlab we shall implement the aformentioned procedure as below. 

```Matlab

parfor r = 1:Rep  % Loop over MC replications (columns first)

   % Define the matrices to store the data
   Tn_matrix = Tn(:,1:end);              

   for j = 1:D   % Loop over designs where D = {1,2,3} 
        %%%%%%%%%%%%%%%%%%%%%% Compute the Statistic %%%%%%%%%%%%%%%%%%%%%%
        % Generate the outcome for the current DGP
        Y = a(j)*Z - b(j)*normpdf(c(j)*Z) + U; 
    
        % Next based on some econometric method, estimate the test statistic Tn 
        %%%%%%%%%%%%%%%%
        % Example below: 
        %%%%%%%%%%%%%%%%
    
        % Coefficient
        Hatbeta  = Hn\Y;   
        Hattheta = X*Hatbeta;     % Unconstrained estimator
        HatU = Y - Hn*Hatbeta;    % Residuals 

        % Compute the test statistic 
        % The calcuation of the test statistic goes here!  
 
        %%%%%%%%%%%%% Obtain Critical Value through Bootstrap %%%%%%%%%%%%%
        BootW = randn(N,Boot);  % Bootstrap weights 

        % The calcuation of the bootsrap test statistic goes here!  
        for i=1:Boot
           
           [BOOTSTRAP PROCEDURE GOES HERE]
           
        end

        % Store the bootstrapped test statistic Tn 
       
           [BOOTSTRAPPED TEST STATISTIC GOES HERE]
       
       %%%%%%%%%%%%%%%%%%%%%%%% Record Rejection %%%%%%%%%%%%%%%%%%%%%%%%%
       Rej(j,r) = ( test_statistic > quantile(BStats,1-alpha));     
    end
end
        
%%%%%%%%%%%%% Compute the Empirical Reject Rates %%%%%%%%%%%%%%%%%
Rej = mean(Rej,2);    
        
```

## 3.3. Econometric Model Fitting and Estimation Examples 

## Example 5



```Matlab






```

# 4. Concluding Remarks

Matlab is a computationally fast programming language which can be employed for various applications in statistics, econometrics, numerical analysis and other. There are plethora of resources for learning programming with Matlab as well as various computational problems worth checking their Matlab implementation or coding your own procedure. 


## References

- Householder, A. S. (1958). Unitary triangularization of a nonsymmetric matrix. Journal of the ACM (JACM), 5(4), 339-342.
- Moravitz Martin, C. D., & Van Loan, C. F. (2007). Solving real linear systems with the complex Schur decomposition. SIAM journal on matrix analysis and applications, 29(1), 177-183.


# Appendix: Parallelism in Matlab


```Matlab

matlabpool('open',8);
tic
start = tic;
clear A
parfor i = 1:100000;
        A(i) = i;
end
stop = toc(start);
stop
matlabpool('close');

```

# Reading List

On Econometric Theory and Applications: 

[1] Davidson, R., & MacKinnon, J. G. (2004). Econometric Theory and Methods (Vol. 5, pp. 189-196). New York: Oxford University Press.

[2] Davidson, J. (2000). Econometric Theory. Blackwell Publishing.


On Matrix Algebra and Numerical Methods: 

[1] Epperson, J. F. (2021). An Introduction to Numerical Methods and Analysis. John Wiley & Sons.

[2] Kharab, A., & Guenther, R. (2018). An Introduction to Numerical Methods: A MATLAB® Approach. CRC press.

[3] Golub, G. H., & Van Loan, C. F. (2013). Matrix computations. JHU press.

## Greek Bibliography:

Γεωργίου, Γ., and Χ. Ξενοφώντος. "Εισαγωγή στο MATLAB." Τμήμα Μαθηματικών και Στατιστικής του Πανεπιστημίου της Κύπρου, Σημειώσεις Εαρινού Εξαμήνου 2007, σελ. 1 306 (2007).


# Disclaimer

The author (Christis G. Katsouris) declares no conflicts of interest.

The proposed Course Syllabus is currently under development and has not been officially undergone quality checks. All rights reserved.

Any errors or omissions are the responsibility of the author.

# Acknowledgements

The author has benefited by participating in workshops and training sessions related to High Performance Computing both at the University of Southampton as well as at University College London (UCL).

Further mathematical/theoretical applications related to Numerical Analysis, Scientific Computing or High Performance Computing can be also found on the website of [Prof. Robert Scheichl](https://katana.iwr.uni-heidelberg.de/people/rob/).

# How to Cite a Website

See: https://www.mendeley.com/guides/web-citation-guide/

# Thank you!
