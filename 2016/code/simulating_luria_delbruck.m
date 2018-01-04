%This script will perform a stochastic simulation of the two predictions for
%the origin of mutation put forth by Max Delbruck and Salvador Luria. 

%We'll start with the case of acquired immunity. We first need to define a few
%parameters. 
n_0 = 100; %Initial number of cells. 
n_gen = 10; %Number of generations. 
mut_rate = 1E-5; %Mutation rate. 
n_sim = 1000; %Number of simulations for good statistics. 

%Under the acquiried immunity hypothesis, cells are only able to 
%mutate when then are exposed to some selective pressure. This makes the
%simulation very easy as we just need to figure out how many cells are present
%at when they are placed under selection and compute the fraction that are
%mutated. We'll loop through each simulation separately. 
n_final = n_0 * 2^n_gen;

for i=1:n_sim
	%Go through each cell and determine if it is mutated. 
%	mutated_cells = rand(n_final, 1)  < mut_rate;
	mutated_cells = 0;
	for j=1:n_final
		%Flip a coin. 
		flip = rand();
		if flip < mut_rate
			mutated_cells = mutated_cells + 1;
		end
	end
	%Now add the number of mutated cells for each simulation to a vector. 
	acquired_immunity(i) = mutated_cells;
end

%It's that easy! Now we can generate a histogram of the mutated cells. 
%Change the bins to range from 0 cells to the maximum number of mutated cells
%taking steps of 1. 
figure(1)
bins = 0:1:max(acquired_immunity);
hist(acquired_immunity,bins)
xlim([0, max(acquired_immunity)]);
xlabel('number of mutated cells')
ylabel('counts')
title('Mutants under Acquired Immunity Hypothesis')


%Now we can consider the random mutation hypothesis. In this case we have to
%keep track of the number of mutated cells at each generation rather than just
%at the end.
for i=1:n_sim
	%Initialize the number of cells at the beginning of each simulation.
	wt_cells = n_0;
	mut_cells = 0;
	%Now go through each generation and find the mutants.
	for j=1:n_gen
		for k = 1:wt_cells;
			%Test if each wild-type cell can mutate.
			flip = rand();
			if flip < mut_rate
				mut_cells = mut_cells + 1;
				%Important to take a wt cell out of the pool.
				wt_cells = wt_cells - 1;
			end
		end

		%Now duplicate the cells. 
		wt_cells = wt_cells * 2;
		mut_cells = mut_cells *2;
	end
	%Now store the number of mutated cells at the end of the simulation. 
	random_mutants(i) = mut_cells;
end

%Now we can show the distribution of mutants for the random mutation
%hypothesis.
figure(2)
bins = 0:1:max(random_mutants);
hist(random_mutants, bins)
xlabel('number of mutated cells')
ylabel('counts')
title('Mutants under Random Mutation Hyptothesis')
xlim([0, max(random_mutants)])


%Now calculate and show the mean and the fano factor for each case. 
acquired_mean = mean(acquired_immunity);
acquired_fano = var(acquired_immunity) / acquired_mean;
random_mean = mean(random_mutants);
random_fano = var(random_mutants) / random_mean;

%Now Display it. 
disp(['Acquired immunity mean is ', num2str(acquired_mean),...
' and fano factor is ', num2str(acquired_fano)])
disp(['Random mutation mean is ', num2str(random_mean),...
' and fano factor is ',num2str(random_fano)])

%The fano factor is huge for the random mutation case, as it should be for a
%non-poissonian process.
