# Introduction-Matlab
Teaching Page (Drafted July 2022).

TIP: Learn, Unlearn & Relearn.

# 1. Preliminary

Matlab is an Object-Oriented Programming Language which has faster computational capabilities in comparison to other Statistical Software. The best way to master a programming language such as Matlab or C++ is to consider reviewing some key concepts from about 50 years ago. Two seminal references are: 

- Tenenbaum, A. M., & Augenstein, M. J. (1986). Data structures using Pascal. Prentice-Hall, Inc.
- Standish, T. A. (1980). Data structure techniques. Addison-Wesley Longman Publishing Co., Inc.

In summary the main idea is to aim to produce clean coding procedures and reproducible Matlab Scripts. This can be helpful when reviewing your code in order to identify bugs. Furthermore, this can be extremely helpful when colloborating on projects with others.  

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%script tutorial.m

format bank 
invest = ('type initial investmentL  ')

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


# 2. Optimization and Estimation Examples

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

```


# 3. Programming in Matlab using functions and programms







# 4. Concluding Remarks


## References

- Householder, A. S. (1958). Unitary triangularization of a nonsymmetric matrix. Journal of the ACM (JACM), 5(4), 339-342.
- Moravitz Martin, C. D., & Van Loan, C. F. (2007). Solving real linear systems with the complex Schur decomposition. SIAM journal on matrix analysis and applications, 29(1), 177-183.


# Appendix: Some useful MATLAB commands







# Reading List

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

Furthermore the author has benefited from the class of [Prof. Robert Scheichl](https://katana.iwr.uni-heidelberg.de/people/rob/) while being an undergraduate student at the Department of Mathematical Sciences of the University of Bath as well as from the class of [Prof. Tassos Magdalinos](https://sites.google.com/site/tmagdalinos/home) while being a PhD student at the Department of Economics of the University of Southampton. 


# How to Cite a Website

See: https://www.mendeley.com/guides/web-citation-guide/
