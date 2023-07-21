# Build 
The example source code has been built on MacOS with Xcode. The program should ok on all Unix-based OS. If you want to run the program, you may put your SMPS files in the spInput folder, and run the  script with
```./LShape_xcode {example_name} ./spInput/ {output_folder}```

If you want to compile the src yourself, a license for IBM cplex is required. You may following this doc to configure IBM CPLEX with Apple Xcode, see https://www.leandro-coelho.com/how-to-configure-ibm-cplex-with-apple-xcode/
And then use the compiled file to run your example.

The src in 'SMIP_VarianceReduction' folder is the algorithm in "Variance Reduced Ensemble Decisions" Section of the paper "Ensemble Variance Reduction Methods for Stochastic Mixed-Integer Programming and their Application to the Stochastic Facility Location Problem".

The src in 'SMIP_OCBA' folder is the algorithm in "Eﬃcient Budget Allocation" Section of the paper "Ensemble Variance Reduction Methods for Stochastic Mixed-Integer Programming and their Application to the Stochastic Facility Location Problem".

If you have any question about the code and data, please contact jiajunx@usc.edu