Written by Rob Hall and Rebecca C. Steorts (2013).
Last updated by Rebecca C. Steorts (2018-04-04). 

1. To run the code, simply download the package.

2. You will need to create a folder called samples/ which will store the output of the MCMC sampler and the matching probabilities as described in Steorts, Hall, Fienberg (2014, 2016). 

3. The main file that runs the code is called "MHSampler_beka2_error_rates_v2.java"

To compile the code, from the console run

javac -classpath mallet.jar:mallet-deps.jar:. MHSampler_beka2_error_rates_v2.java
To run the code and output to a log file 

java -classpath mallet.jar:mallet-deps.jar:. MHSampler_beka2_error_rates_v2 > logfile.txt

Remark: MHSampler_beka2_matches_v2.java allows you to change the data sets and settings of the Gibbs sampler (number of Gibbs iterations, burn-in, etc). 

4. The rest of the output can be processed using R files or other scripts. 

Specifically, you need a script in order to put all the MCMC samples together to form the linkage structure as talked about in the paper. Once this has been completed, one can look at diagnostic plots and posterior analysis. 

Using the matching file outputted from this algorithm, one can look at the MPMMS and posterior matching probabilities. 

5. The perl script "consolidate_samples_uniform_filenames.pl" puts together the Gibbs samples "linkage_iter" and writes them to a file called "posterior.link"